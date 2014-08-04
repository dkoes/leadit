/*
 * DatabaseCreator.cpp
 *
 *  Created on: Jul 16, 2014
 *      Author: dkoes
 */

#include <DatabaseCreator.h>
#include <fstream>
#include <rdkit/GraphMol/FileParsers/MolSupplier.h>
#include <rdkit/GraphMol/FileParsers/MolWriters.h>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>

//creation
DatabaseCreator::DatabaseCreator(const vector<filesystem::path>& dbs,
		const Reaction& react, const DatabaseConfiguration& config) : dbpaths(dbs),
		rxn(react), valid(false)
{
	scaffoldIndex.initialize(rxn, config.scaffoldRMSDcutoff, config.connectPointCutoff);

	//make the directories
	for(unsigned i = 0, n = dbpaths.size(); i < n; i++)
	{
		if(!filesystem::create_directory(dbpaths[i]))
			return;
	}
	valid = true;
}

//addition
DatabaseCreator::DatabaseCreator(const vector<filesystem::path>& dbs): dbpaths(dbs)
{
	assert(dbpaths.size() > 0);

	filesystem::path rxnpath = dbpaths[0] / "rxninfo";
	ifstream rxnf(rxnpath.c_str());
	if(!rxnf) return;
	rxn.read(rxnf);

	filesystem::path idxpath = dbpaths[0] / "scaffoldIdx";
	ifstream scaffin(idxpath.c_str());
	if(!scaffin) return;
	scaffoldIndex.read(scaffin);

	unsigned sz = scaffoldIndex.size();
	fragments.resize(sz);
	abort(); //need to read in frags
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
	for(unsigned s = 0, ns = fragments.size(); s < ns; s++)
	{
		unsigned long possibleconfs = 1;
		for(unsigned p = 0, np = fragments[s].size(); p < np; p++)
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
	for(unsigned s = 0, ns = fragments.size(); s < ns; s++)
	{
		unsigned possibleconfs = 1;
		unsigned possiblemols = 1;
		out << s << "\t" << scaffoldIndex.getCnt(s) << "\t";
		for(unsigned p = 0, np = fragments[s].size(); p < np; p++)
		{
			const FragmentIndexer& fi = fragments[s][p];
			out << fi.numFragments() << " (" << fi.numFragmentConformers() << ")\t";
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
	if(molfile.extension() == ".gz")
	{
		in.push(iostreams::gzip_decompressor());
	}
	in.push(inmols);

	MCForwardSDMolSupplier supplier(&in, false);
	vector<MOL_SPTR_VECT> pieces;
	vector<ROMOL_SPTR> core;

	int cnt = 0;
	while(!supplier.atEnd())
	{
		ROMOL_SPTR mol = ROMOL_SPTR(supplier.next());
		rxn.decompose(*mol, pieces, core);
		assert(pieces.size() == core.size());

		//treat each match separately - highly similar confs will get weeded out anyway
		for(unsigned i = 0, n = core.size(); i < n; i++)
		{
			ROMOL_SPTR coremol = core[i];
			for(unsigned c = 0, nc = coremol->getNumConformers(); c < nc; c++)
			{
				Conformer& conf = coremol->getConformer(c);
				//categorize the scaffold
				Orienter orient;
				unsigned sindex = scaffoldIndex.addScaffold(conf, orient);
				//position conformer to be aligned to core scaffold
				orient.reorient(conf.getPositions());

				//check for new scaffold
				if(fragments.size() == sindex)
					fragments.push_back(vector<FragmentIndexer>(pieces[i].size(), FragmentIndexer(config.reactantRMSDcutoff)));
				assert(sindex < fragments.size());
				//each scaffold_reactantpos is a unique database
				for(unsigned p = 0, np = pieces[i].size(); p < np; p++)
				{
					ROMOL_SPTR frag = pieces[i][p];
					assert(frag->getNumConformers() == nc);
					Conformer& fragconf = frag->getConformer(c);
					orient.reorient(fragconf.getPositions()); //align to match scaffold

					fragments[sindex][p].add(fragconf);
				}

				if(verbose)
				{
					cnt++;
					cout << "CNTS " << cnt << "\t" << totalConformers() << "\t" << fragments.size() << "\n";
				}
			}
		}
	}
	cout << "Index size " << scaffoldIndex.size() << "\n";
	dumpCounts(cout);
}


//generate indices and write all data to disk
void DatabaseCreator::finalize()
{

}
