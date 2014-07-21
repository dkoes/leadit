/*
 * 5/21/2014
 * David Koes (dkoes@pitt.edu)
 *
 * The intent for leadmaker is that it will enable lead optimization of existing
 * hits by deconstructing a ligand's chemistry and, using the same reaction pathway,
 * reconstruct alternative leads that match specified shape and pharmacophore criteria.
 *
 * leadmaker constructs databases of molecular fragments.  Each database is
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

using namespace std;

enum CommandEnum
{
	CreateDatabase,
	AddMolecules,
	SearchDatabase,
	DatabaseInfo,
	Server
};

cl::opt<CommandEnum> Command(cl::desc("Operation to perform:"), cl::Required,
		cl::values(clEnumVal(CreateDatabase, "Create a new reaction database"),
				clEnumVal(AddMolecules,
						"Add conformers to database (regenerates indices)"),
				clEnumVal(SearchDatabase, "Search database for leadmaker query"),
				clEnumVal(DatabaseInfo, "Print database information"),
				clEnumVal(Server, "Start leadmaker server"),
				clEnumValEnd));
cl::list<string> Databases("dbdir", cl::desc("database directory(s)"));
cl::list<string> inputFiles("in", cl::desc("input file(s)"));
cl::list<string> outputFiles("out", cl::desc("output file(s)"));

cl::opt<string> reactionFile("rxn", cl::desc("reaction SMARTS file"));


//create a database using command line arguments
static void handle_create()
{
	if(Databases.size() == 0)
	{
		cerr << "Require database for create\n";
		exit(-1);
	}

	vector<filesystem::path> dbpaths;
	for(unsigned i = 0, n = Databases.size(); i < n; i++)
		dbpaths.push_back(filesystem::path(Databases[i]));

	//reaction file
	filesystem::path rxnf(reactionFile);
	if(!filesystem::exists(rxnf))
	{
		cerr << rxnf << " reaction file does not exist\n";
		exit(-1);
	}
	Reaction rxn(rxnf);
	if(!rxn.isValid())
	{
		cerr << "Invalid reaction\n";
		exit(-1);
	}
	cout << rxn;
	//open database for creation
	DatabaseCreator dbcreator(dbpaths, rxn);

	if(!dbcreator.isValid())
	{
		cerr << "Error creating database\n";
		exit(-1);
	}
	for(unsigned i = 0, n = inputFiles.size(); i < n; i++)
	{
		filesystem::path infile(inputFiles[i]);
		if(!filesystem::exists(infile))
		{
			cerr << infile << " does not exists. Skipping.\n";
			continue;
		}
		dbcreator.add(infile);
	}

	//create indices
	dbcreator.finalize();
}

//append to database using command line argument values
static void handle_add()
{
	if(Databases.size() == 0)
	{
		cerr << "Require database for add\n";
		exit(-1);
	}

	vector<filesystem::path> dbpaths;
	for(unsigned i = 0, n = Databases.size(); i < n; i++) {
		filesystem::path dbpath = Databases[i];
		if(!filesystem::exists(dbpath))
		{
			cerr << dbpath << " does not exist\n";
			exit(-1);
		}
		dbpaths.push_back(dbpath);
	}

	//open database for appending
	DatabaseCreator dbcreator(dbpaths);


	if(!dbcreator.isValid())
	{
		cerr << "Error opening database\n";
		exit(-1);
	}

	for(unsigned i = 0, n = inputFiles.size(); i < n; i++)
	{
		filesystem::path infile(inputFiles[i]);
		if(!filesystem::exists(infile))
		{
			cerr << infile << " does not exists. Skipping.\n";
			continue;
		}
		dbcreator.add(infile);
	}

	//create indices
	dbcreator.finalize();
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
		case SearchDatabase: //f
		case DatabaseInfo: //f
		case Server: //f
		default:
			cerr << "Command not yet implemented\n";
			return -1;
	}

	return 0;
}
