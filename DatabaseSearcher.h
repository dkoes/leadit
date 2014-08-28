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
#include "QueryObject.h"

using namespace boost;

class DatabaseSearcher
{
	vector<boost::filesystem::path> dbpaths;
	Reaction rxn;

	ScaffoldIndexer scaffoldIndex;

	vector<vector<vector<FragmentSearcher> > > fragments; //index first by scaffold position, then by reactant position, then by dir
	bool valid;

public:

	//fully specify the location of a result, trade redundancy for simplicity
	struct Result
	{
		unsigned scaffoldPos;
		unsigned reactantPos;
		unsigned dir;
		FragmentSearcher::Result fragRes;

		Result() :
				scaffoldPos(0), reactantPos(0), dir(0)
		{
		}

		Result(unsigned s, unsigned r, unsigned d,
				const FragmentSearcher::Result& fr) :
				scaffoldPos(s), reactantPos(r), dir(d), fragRes(fr)
		{
		}
	};

	DatabaseSearcher(const vector<filesystem::path>& dbpaths);
	~DatabaseSearcher()
	{
	}

	//true if no problems setting up creation
	bool isValid() const
	{
		return valid;
	}

	const Reaction& getReaction() const
	{
		return rxn;
	}
	unsigned long totalConformers() const;

	void search(ROMOL_SPTR ref, unsigned reactant, const QueryObject& small,
			const QueryObject& big, vector<Result>& results); //TODO pharma

	//write result sdf to out
	void writeSDF(const Result& res, ostream& out) const;

};

#endif /* DATABASESEARCHER_H_ */
