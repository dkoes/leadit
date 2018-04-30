/*
 * 5/21/2014
 * David Koes (dkoes@pitt.edu)
 *
 * The intent for leadit is that it will enable lead optimization of existing
 * hits by deconstructing a ligand's chemistry and, using the same reaction pathway,
 * reconstruct alternative leads that match specified shape and pharmacophore criteria.
 *
 * leadit constructs databases of molecular fragments.  Each database is
 * constructed around a single reaction framework.  It is populated using
 * conformers that are compatible with that framework.  A minimal scaffold is
 * identified from the reaction.  Scaffold conformers are clustered with respect
 * to RMSD and connection point deviations.  Reaction components are clustered
 * with respect to unminimized RMSD after alignment to the appropriate scaffold.
 * Only representative fragments are retained.
 *
 * The fragment conformations will be indexed using shapedb with some color
 * information.  I may also add in a kd-tree lookup if the color shapedb color
 * scheme doesn't work as well as I hope.
 *
 * It will also be possible to print statistics of the current database
 * and to deconstruct single compounds.
 *
 * Search functionality for single and multiple databases will be implemented.
 *
 * I will need to add robust support for I/O parallelism.
 *
 * Eventually a server mode will be added.
 */

#include "CommandLine2/CommandLine.h"
#include <string>
#include <iostream>

#include "Reaction.h"
#include "DatabaseCreator.h"
#include "DatabaseSearcher.h"
#include <GraphMol/FileParsers/MolWriters.h>

using namespace std;

enum CommandEnum
{
	CreateDatabase,
	AddMolecules,
	SearchDatabase,
	DatabaseInfo,
	LigandInfo,
	Server
};

cl::opt<CommandEnum> Command(cl::desc("Operation to perform:"), cl::Required,
		cl::values(clEnumVal(CreateDatabase, "Create a new reaction database"),
				clEnumVal(AddMolecules,
						"Add conformers to database (regenerates indices)"),
				clEnumVal(SearchDatabase,
						"Search database for leadit query"),
				clEnumVal(DatabaseInfo, "Print database information"),
				clEnumVal(LigandInfo,
						"Print decomposition of passed ligand(s)"),
				clEnumVal(Server, "Start leadit server"),
				clEnumValEnd));
cl::list<string> Databases("dbdir", cl::desc("database directory(s)"));
cl::list<string> inputFiles("in", cl::desc("input file(s)"));
cl::opt<string> outputFile("out", cl::desc("output file"));

cl::opt<string> IncludeMol("inc-mol",
		cl::desc("Ligand minimum volume shape constraint."));
cl::opt<double> IncludeShrink("inc-shrink",
		cl::desc("Amount to reduce ligand minimum shape."), cl::init(0));
cl::opt<string> ExcludeMol("exc-mol",
		cl::desc("Receptor excluded shape constraint."));
cl::opt<double> ExcludeShrink("exc-shrink",
		cl::desc("Amount to reduce excluded shape."), cl::init(0));
cl::opt<string> PharmaQuery("pharma",
		cl::desc("Pharmacophore search constraints"));

cl::opt<int> ReactantPos("rpos",
		cl::desc("Reactant position to replace for search"), cl::init(-1));
cl::opt<string> RefLigand("ref", cl::desc("Reference scaffold for search"));

cl::opt<string> reactionFile("rxn", cl::desc("reaction SMARTS file"));

cl::opt<bool> Force("force", cl::desc("Overwrite any existing database"),
		cl::init(false));
cl::opt<double> ScaffoldRMSD("scaffold-rmsd",
		cl::desc("Maximum RMSD for scaffolds to be merged"), cl::init(0.5));
cl::opt<double> ConnectCutoff("connect-cutoff",
		cl::desc("Maximum distance allowed between connection points"),
		cl::init(0.1));
cl::opt<double> ReactantRMSD("reactant-rmsd",
		cl::desc("Maximum RMSD for reactants to be merged"), cl::init(0.5));

cl::opt<bool> Verbose("verbose", cl::desc("verbose output"));

//create a database using command line arguments
static void handle_create()
{
	if (Databases.size() == 0)
	{
		cerr << "Require database for create\n";
		exit(-1);
	}

	vector<filesystem::path> dbpaths;
	for (unsigned i = 0, n = Databases.size(); i < n; i++)
		dbpaths.push_back(filesystem::path(Databases[i]));

	//reaction file
	filesystem::path rxnf(reactionFile);
	if (!filesystem::exists(rxnf))
	{
		cerr << rxnf << " reaction file does not exist\n";
		exit(-1);
	}
	Reaction rxn(rxnf);
	if (!rxn.isValid())
	{
		cerr << "Invalid reaction\n";
		exit(-1);
	}
	cout << rxn;
	//parameters
	DatabaseCreator::DatabaseConfiguration config;
	config.force = Force;
	config.reactantRMSDcutoff = ReactantRMSD;
	config.connectPointCutoff = ConnectCutoff;
	config.scaffoldRMSDcutoff = ScaffoldRMSD;

	//open database for creation
	DatabaseCreator dbcreator(dbpaths, rxn, config);

	if (!dbcreator.isValid())
	{
		cerr << "Error creating database\n";
		exit(-1);
	}
	for (unsigned i = 0, n = inputFiles.size(); i < n; i++)
	{
		filesystem::path infile(inputFiles[i]);
		if (!filesystem::exists(infile))
		{
			cerr << infile << " does not exists. Skipping.\n";
			continue;
		}
		dbcreator.add(infile, Verbose);
	}

	//create indices
	dbcreator.finalize();
}

//append to database using command line argument values
static void handle_add()
{
	if (Databases.size() == 0)
	{
		cerr << "Require database for add\n";
		exit(-1);
	}

	vector<filesystem::path> dbpaths;
	for (unsigned i = 0, n = Databases.size(); i < n; i++)
	{
		filesystem::path dbpath = Databases[i];
		if (!filesystem::exists(dbpath))
		{
			cerr << dbpath << " does not exist\n";
			exit(-1);
		}
		dbpaths.push_back(dbpath);
	}

	//open database for appending
	DatabaseCreator dbcreator(dbpaths);

	if (!dbcreator.isValid())
	{
		cerr << "Error opening database\n";
		exit(-1);
	}

	for (unsigned i = 0, n = inputFiles.size(); i < n; i++)
	{
		filesystem::path infile(inputFiles[i]);
		if (!filesystem::exists(infile))
		{
			cerr << infile << " does not exists. Skipping.\n";
			continue;
		}
		dbcreator.add(infile);
	}

	//create indices
	dbcreator.finalize();
}

static void handle_dbinfo()
{
	if (Databases.size() == 0)
	{
		cerr << "Require database for dbinfo\n";
		exit(-1);
	}

	vector<filesystem::path> dbpaths;
	for (unsigned i = 0, n = Databases.size(); i < n; i++)
	{
		filesystem::path dbpath = Databases[i];
		if (!filesystem::exists(dbpath))
		{
			cerr << dbpath << " does not exist\n";
			exit(-1);
		}
		dbpaths.push_back(dbpath);
	}

	//open database for appending
	DatabaseSearcher dbsearcher(dbpaths);

	if (!dbsearcher.isValid())
	{
		cerr << "Error opening database\n";
		exit(-1);
	}

	cout << dbsearcher.getReaction() << "\n";
	cout << dbsearcher.totalConformers() << " total conformers available\n";
}

//print out information on all ligands in ligfile
//if out is valid,
static void print_ligand_info(Reaction& rxn, const string& ligfile,
		ostream& out)
{
	ifstream inmols(ligfile.c_str());
	iostreams::filtering_stream<iostreams::input> in;
	if (filesystem::extension(ligfile) == ".gz")
	{
		in.push(iostreams::gzip_decompressor());
	}
	in.push(inmols);
	MCForwardSDMolSupplier supplier(&in, false);
	vector<MOL_SPTR_VECT> pieces;
	vector<ROMOL_SPTR> core;

	SDWriter writer(&out);
	while (!supplier.atEnd())
	{
		ROMOL_SPTR mol = ROMOL_SPTR(supplier.next());
		ROMOL_SPTR m(MolOps::addHs(*mol)); //for proper match must have hydrogens
		rxn.decompose(m, pieces, core);
		assert(core.size() == pieces.size());
		for (unsigned i = 0, n = core.size(); i < n; i++)
		{
			cout << MolToSmiles(*core[i]) << "\t";
			if (out)
			{
				//hack to quickly set name
				out << "core ";
				writer.write(*core[i]);
			}
			for (unsigned j = 0, m = pieces[i].size(); j < m; j++)
			{
				cout << MolToSmiles(*pieces[i][j]) << " ";
				if (out)
				{
					//quickly set name
					out << "piece" << j << " ";
					writer.write(*pieces[i][j]);
				}
			}
			cout << "\n";
		}
	}
}

//perform ligand decomposition
static void handle_ligand_info()
{
	if (inputFiles.size() == 0)
	{
		cerr << "Need ligand input\n";
		exit(-1);
	}

	ofstream out(outputFile.c_str());

	if (Databases.size() > 0)
	{
		//if databases are specified, look for a matching reaction
		BOOST_FOREACH(const string& dbfile, Databases)
		{
			filesystem::path rxnf = filesystem::path(dbfile) / "rxninfo";
			ifstream rxnin(rxnf.c_str());
			Reaction rxn;
			rxn.read(rxnin);
			BOOST_FOREACH(const string& ligfile, inputFiles)
			{
				print_ligand_info(rxn, ligfile, out);
			}
		}
	}
	else if (reactionFile.size() > 0)
	{
		filesystem::path rxnf(reactionFile);
		ifstream rxnin(rxnf.c_str());
		Reaction rxn;
		rxn.read(rxnin);
		BOOST_FOREACH(const string& ligfile, inputFiles)
		{
			print_ligand_info(rxn, ligfile, out);
		}
	}
	else
	{
		cerr << "Need database or reaction file to analyze ligand\n";
		exit(-1);
	}
}

//read a single molecule from the file (assumed sdf)
static ROMOL_SPTR readOneMol(const string& filename)
{
	ifstream inmols(filename.c_str());
	iostreams::filtering_stream<iostreams::input> in;
	if (filesystem::extension(filename) == ".gz")
	{
		in.push(iostreams::gzip_decompressor());

		string stripped = filename.substr(0,filename.find_last_of("."));
		if (filesystem::extension(stripped) != ".sdf")
		{
	    cerr << "Sorry, currently only sdf files are supported, " << filename << " does not appear to be an sdf.\n";
	    exit(-1);
		}
	}
	else if (filesystem::extension(filename) != ".sdf")
	{
	  cerr << "Sorry, currently only sdf files are supported, " << filename << " does not appear to be an sdf.\n";
	  exit(-1);
	}
	in.push(inmols);

	ForwardSDMolSupplier molreader(&in, false);
	ROMOL_SPTR mol(molreader.next());

	return mol;
}

static void handle_search()
{
	if (ReactantPos < 0)
	{
		cerr << "Need position of reactant to replace\n";
		exit(-1);
	}

	if (RefLigand.size() == 0)
	{
		cerr << "Need reference ligand for scaffold positioning\n";
		exit(-1);
	}

	if (Databases.size() == 0)
	{
		cerr << "Need databases to search\n";
		exit(-1);
	}

	//setup databases - must all be the same, just striped
	vector<filesystem::path> dbpaths; dbpaths.reserve(Databases.size());
	BOOST_FOREACH(const string& dbfile, Databases)
	{
		dbpaths.push_back(filesystem::path(dbfile));
	}

	DatabaseSearcher searcher(dbpaths);

	if(!searcher.isValid())
	{
		cerr << "Problem reading in database\n";
		exit(-1);
	}

	//read in reference ligand
	ROMOL_SPTR refmol = readOneMol(RefLigand);

	MolecularQueryObject small, big(true);

	if(IncludeMol.size() > 0)
	{
		ROMOL_SPTR smalllig = readOneMol(IncludeMol);
		small = MolecularQueryObject(*smalllig, IncludeShrink, false);
	}
	if(ExcludeMol.size() > 0)
	{
		ROMOL_SPTR biglig = readOneMol(ExcludeMol);
		big = MolecularQueryObject(*biglig, ExcludeShrink, true);
	}


	vector<DatabaseSearcher::Result> results;
	searcher.search(refmol, ReactantPos, small, big, results);

	cout << results.size() << " hits\n";

	if(outputFile.size() > 0)
	{
		ofstream out(outputFile.c_str());
		BOOST_FOREACH(DatabaseSearcher::Result& r, results)
		{
			searcher.writeSDF(r, out);
		}
	}
}

int main(int argc, char *argv[])
{
	cl::ParseCommandLineOptions(argc, argv);

	switch (Command)
	{
		case CreateDatabase:
			handle_create();
			break;
		case AddMolecules:
			handle_add();
			break;
		case DatabaseInfo:
			handle_dbinfo();
			break;
		case LigandInfo:
			handle_ligand_info();
			break;
		case SearchDatabase:
			handle_search();
			break;
		case Server: //f
		default:
			cerr << "Command not yet implemented\n";
			return -1;
	}

	return 0;
}
