/*
 * DatabaseSearcher.h
 *
 *  Class for reading in a leadmaker database for searching and ligand analysis.
 *  Created on: Aug 26, 2014
 *      Author: dkoes
 */

#ifndef DATABASESEARCHER_H_
#define DATABASESEARCHER_H_

#include <boost/filesystem.hpp>

#include "Reaction.h"
#include "ScaffoldIndexer.h"
#include "FragmentSearcher.h"

using namespace boost;

class DatabaseSearcher
{
	vector<boost::filesystem::path> dbpaths;
	Reaction rxn;

	ScaffoldIndexer scaffoldIndex;

	vector< vector< vector<FragmentSearcher> > > fragments; //index first by scaffold position, then by reactant position, then by dir
	bool valid;

public:
	DatabaseSearcher(const vector<filesystem::path>& dbpaths);
	~DatabaseSearcher() {}

	//true if no problems setting up creation
	bool isValid() const { return valid; }

	const Reaction& getReaction() const { return rxn; }
	unsigned long totalConformers() const;

};

#endif /* DATABASESEARCHER_H_ */
