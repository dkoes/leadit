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

//add conformers in molfile, only adds data, does not generate indices
void DatabaseCreator::add(const filesystem::path& molfile)
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

	while(!supplier.atEnd())
	{
		ROMOL_SPTR mol = ROMOL_SPTR(supplier.next());
		rxn.decompose(*mol, pieces, core);
		assert(pieces.size() == core.size());
		//treat each match separately - highly similar confs will get weeded out anyway
		for(unsigned i = 0, n = core.size(); i < n; i++)
		{
			ROMOL_SPTR coremol = core[i];
			for(ROMol::ConstConformerIterator itr = coremol->beginConformers(), end = coremol->endConformers();
					itr != end; ++itr)
			{
				const CONFORMER_SPTR conf = *itr;
				//categorize the scaffold
				Orienter orient;
				unsigned sindex = scaffoldIndex.addScaffold(conf, orient);
				//position conformer to be aligned to core scaffold
				orient.reorient(conf->getPositions());

				//each scaffold_reactantpos is a unique database
			}
		}
	}
	cout << "Index size " << scaffoldIndex.size() << "\n";
	scaffoldIndex.dumpCounts(cout);
}


//generate indices and write all data to disk
void DatabaseCreator::finalize()
{

}
