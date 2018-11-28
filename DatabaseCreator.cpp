/*
 * DatabaseCreator.cpp
 *
 *  Created on: Jul 16, 2014
 *      Author: dkoes
 */

#include "DatabaseCreator.h"
#include <fstream>
#include <rdkit/GraphMol/FileParsers/MolSupplier.h>
#include <rdkit/GraphMol/FileParsers/MolWriters.h>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include "pharmacophores.h"
//creation
DatabaseCreator::DatabaseCreator(const vector<filesystem::path>& dbs,
		const Reaction& react, const DatabaseConfiguration& config) :
		dbpaths(dbs),
				rxn(react), valid(false)
{
	scaffoldIndex.initialize(rxn, config.scaffoldRMSDcutoff,
			config.connectPointCutoff);

	//make the directories
	for (unsigned i = 0, n = dbpaths.size(); i < n; i++)
	{
		if(filesystem::exists(dbpaths[i]) && config.force)
		{
			//try to remove any existing directories
			filesystem::remove_all(dbpaths[i]);
		}
		if (!filesystem::create_directory(dbpaths[i]))
			return;
	}
	valid = true;
}

//addition
DatabaseCreator::DatabaseCreator(const vector<filesystem::path>& dbs) :
		dbpaths(dbs)
{
	assert(dbpaths.size() > 0);

	filesystem::path rxnpath = dbpaths[0] / "rxninfo";
	ifstream rxnf(rxnpath.c_str());
	if (!rxnf)
	{
		cerr << "Could not read " << rxnpath << "\n";
		return;
	}
	rxn.read(rxnf);

	filesystem::path idxpath = dbpaths[0] / "scaffoldIdx";
	ifstream scaffin(idxpath.c_str());
	if (!scaffin)
	{
		cerr << "Could not read " << idxpath << "\n";
		return;
	}
	scaffoldIndex.read(scaffin);

	unsigned sz = scaffoldIndex.size();
	unsigned nr = rxn.numPieces();
	fragments.resize(sz);

	for(unsigned f = 0; f < sz; f++)
	{
		fragments[f].resize(nr);
		for(unsigned p = 0; p < nr; p++)
		{
			vector<filesystem::path> fragdirs;
			for(unsigned i = 0, n = dbpaths.size(); i < n; i++)
			{
				filesystem::path fragdir = dbpaths[i] / lexical_cast<string>(f) / lexical_cast<string>(p);
				if(!filesystem::is_directory(fragdir))
				{
					cerr << "Could not read " << fragdir << "\n";
					return;
				}
				fragdirs.push_back(fragdir);
			}
			bool success = fragments[f][p].read(fragdirs);
			if(!success)
				return;
		}
	}

	valid = true;
}

DatabaseCreator::~DatabaseCreator()
{
}

//true if no problems setting up creation
bool DatabaseCreator::isValid() const
{
	return valid;
}

unsigned long DatabaseCreator::totalConformers() const
{
	unsigned long totalconfs = 0;
	for (unsigned s = 0, ns = fragments.size(); s < ns; s++)
	{
		unsigned long possibleconfs = 1;
		for (unsigned p = 0, np = fragments[s].size(); p < np; p++)
		{
			const FragmentIndexer& fi = fragments[s][p];
			possibleconfs *= fi.numFragmentConformers();
		}

		totalconfs += possibleconfs;
	}
	return totalconfs;
}

void DatabaseCreator::dumpCounts(ostream& out) const
		{
	unsigned long totalconfs = 0;
	for (unsigned s = 0, ns = fragments.size(); s < ns; s++)
	{
		unsigned possibleconfs = 1;
		unsigned possiblemols = 1;
		out << s << "\t" << scaffoldIndex.getCnt(s) << "\t";
		for (unsigned p = 0, np = fragments[s].size(); p < np; p++)
		{
			const FragmentIndexer& fi = fragments[s][p];
			out << fi.numFragments() << " (" << fi.numFragmentConformers()
					<< ")\t";
			possibleconfs *= fi.numFragmentConformers();
			possiblemols *= fi.numFragments();
		}

		out << possiblemols << " " << possibleconfs;
		totalconfs += possibleconfs;
		out << "\n";
	}
	cout << "Total possible structures: " << totalconfs << "\n";
}

//add conformers in molfile, only adds data, does not generate indices
void DatabaseCreator::add(const filesystem::path& molfile, bool verbose /* = false */)
{
	//handle gzipped sdf
	ifstream inmols(molfile.c_str());
	iostreams::filtering_stream<iostreams::input> in;
	if (molfile.extension() == ".gz")
	{
		in.push(iostreams::gzip_decompressor());
	}
	in.push(inmols);

	MCForwardSDMolSupplier supplier(&in, false);

	int cnt = 0;
	while (!supplier.atEnd())
	{
		ROMOL_SPTR mol;
		ROMOL_SPTR m;
		try {
			mol = ROMOL_SPTR(supplier.next());
			if(!mol) throw std::runtime_error("Couldn't parse mol");
			m = ROMOL_SPTR(MolOps::addHs(*mol,false,true)); //for proper match must have hydrogens
		} catch (const std::exception& e) {
			std::cerr << "WARNING: Could not process molecule ";
			if(mol) {
				std::cerr << RDKit::MolToSmiles(*mol) << "\n";
			} else {
				std::cerr << "\n";
			}
			continue;
		}
		//compute pharmacophore features of mol and annotate atoms with feature types
		assignPharmacophoreAtomProperties(m);
		vector<Reaction::Decomposition> decomps;
		if(!rxn.decompose(m, decomps))
		{
		  std::cerr << "WARNING: Could not decompose " << RDKit::MolToSmiles(*m) << "\n";
		  continue;
		}

		//treat each match separately - highly similar confs will get weeded out anyway
		for (unsigned i = 0, n = decomps.size(); i < n; i++)
		{
			const Reaction::Decomposition& d = decomps[i];
			//extract the fragments
			MOL_SPTR_VECT frags;
			vector< vector<FragBondInfo> > fragbonds;
			fragbonds.resize(d.pieces.size());
			bool failed = false;
			for (unsigned p = 0, np = d.pieces.size(); p < np; p++)
			{
				ROMOL_SPTR frag = d.extractPiece(m, p, fragbonds[p]); //labels connecting atoms
				if(!frag) {
					failed = true;
					break;
				}
				frags.push_back(frag);
			}
			//be robust to errors, even though they shouldn't happen
			if(failed) {
				std::cerr << "ERROR: Failed to process decomposition of " << RDKit::MolToSmiles(*m) << "\n";
				continue;
			}

			for (unsigned c = 0, nc = m->getNumConformers(); c < nc; c++)
			{
				Conformer& conf = m->getConformer(c);
				//categorize the scaffold
				Orienter orient;
				unsigned sindex = scaffoldIndex.addScaffold(conf, d, orient);
				//position conformer to be aligned to core scaffold
				orient.reorient(conf.getPositions()); //this probably isn't necessary, fragments get reoriented themselves

				//check for new scaffold
				if (fragments.size() == sindex)
					fragments.push_back(vector<FragmentIndexer>(d.pieces.size(),
									FragmentIndexer(config.reactantRMSDcutoff)));
				assert(sindex < fragments.size());
				//each scaffold_reactantpos is a unique database
				for (unsigned p = 0, np = d.pieces.size(); p < np; p++)
				{
					ROMOL_SPTR frag = frags[p];
					assert(frag->getNumConformers() == nc);
					Conformer& fragconf = frag->getConformer(c);
					orient.reorient(fragconf.getPositions()); //align to match scaffold

					fragments[sindex][p].add(fragconf, fragbonds[p]);
				}

				if (verbose)
				{
					cnt++;
					cout << "CNTS " << cnt << "\t" << totalConformers() << "\t"
							<< fragments.size() << "\n";
				}
			}
		}
	}
	cout << "Index size " << scaffoldIndex.size() << "\n";
	dumpCounts(cout);
}

//generate indices and write all data to disk
//distributes data across database directories
void DatabaseCreator::finalize()
{
	for (unsigned i = 0, n = dbpaths.size(); i < n; i++)
	{
		//duplicate rxinfo and scaffold idx
		filesystem::path rxnpath = dbpaths[i] / "rxninfo";
		ofstream rxnf(rxnpath.c_str());
		if (!rxnf)
		{
			cerr << "Could not write " << rxnpath << "\n";
			abort();
		}
		rxn.write(rxnf);

		filesystem::path idxpath = dbpaths[i] / "scaffoldIdx";
		ofstream scaffout(idxpath.c_str());
		if (!scaffout)
		{
			cerr << "Could not write " << idxpath << "\n";
			abort();
		}
		scaffoldIndex.write(scaffout);
	}

	//the molecule data in fragments gets striped
	//first setup the directories
	vector<vector<vector<filesystem::path> > > fragdirs(fragments.size());
	for (unsigned f = 0, nf = fragments.size(); f < nf; f++)
	{
		fragdirs[f].resize(fragments[f].size());

		//create directory corresponding to scaffold index if necessary
		for (unsigned d = 0, nd = dbpaths.size(); d < nd; d++)
		{
			filesystem::path fragidxdir = dbpaths[d] / lexical_cast<string>(f);
			if (!filesystem::is_directory(fragidxdir)
					&& !filesystem::create_directory(fragidxdir))
			{
				cerr << "Could not make directory " << fragidxdir << "\n";
				abort();
			}
			//make directory for each reactant pos
			for (unsigned p = 0, np = fragments[f].size(); p < np; p++)
			{
				filesystem::path posdir = fragidxdir / lexical_cast<string>(p);
				if (!filesystem::is_directory(posdir)
						&& !filesystem::create_directory(posdir))
				{
					cerr << "Could not make directory " << fragidxdir << "\n";
					abort();
				}
				fragdirs[f][p].push_back(posdir);
			}
		}
	}

	for (unsigned f = 0, nf = fragments.size(); f < nf; f++)
	{
		for (unsigned p = 0, np = fragments[f].size(); p < np; p++)
		{
			fragments[f][p].write(fragdirs[f][p]);
		}
	}
}
