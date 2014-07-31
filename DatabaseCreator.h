/*
 * DatabaseCreator.h
 *
 *  Created on: Jul 16, 2014
 *      Author: dkoes
 *
 *  This class contains the data structures and algorithms necessary
 *  to create a searchable leadmaker database.
 */

#ifndef DATABASECREATOR_H_
#define DATABASECREATOR_H_

#include <boost/filesystem.hpp>

#include "Reaction.h"
#include "ScaffoldIndexer.h"
#include "FragmentIndexer.h"

using namespace boost;

class DatabaseCreator
{
public: //class types
	//various parameters for clustering and so on
	struct DatabaseConfiguration
	{
		double scaffoldRMSDcutoff; //maximum rmsd for scaffolds to be identical
		double connectPointCutoff; //maximum allowed deviation of connecting points
		double reactantRMSDcutoff; //maximum rmsd for reactant conformation to be considered novel
		//sensible defaults
		DatabaseConfiguration() :
				scaffoldRMSDcutoff(0.5), connectPointCutoff(0.1),
						reactantRMSDcutoff(0.5)
		{
		}
	};

private:
	vector<boost::filesystem::path> dbpaths;
	Reaction rxn;
	DatabaseConfiguration config;

	ScaffoldIndexer scaffoldIndex;

	vector< vector<FragmentIndexer> > fragments; //index first by scaffold position, then by reactant position
	bool valid;

public:
	//create a database for a specific reaction
	DatabaseCreator(const vector<filesystem::path>& dbpaths, const Reaction& react,
			const DatabaseConfiguration& config = DatabaseConfiguration());
	//open an existing database for appending
	DatabaseCreator(const vector<filesystem::path>& dbpaths);

	~DatabaseCreator();

	//true if no problems setting up creation
	bool isValid() const;

	//add conformers in molfile, only adds data, does not generate indices
	void add(const filesystem::path& molfile);

	//generate indices and write all data to disk
	void finalize();

	void dumpCounts(ostream& out) const;
};

#endif /* DATABASECREATOR_H_ */
