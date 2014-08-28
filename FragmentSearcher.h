/*
 * FragmentSearcher.h
 *
 *  Created on: Aug 26, 2014
 *      Author: dkoes
 */

#ifndef FRAGMENTSEARCHER_H_
#define FRAGMENTSEARCHER_H_

#include <boost/filesystem.hpp>
#include <vector>
#include <string>
#include "MolMatcher.h"
#include "shapedb/GSSTreeSearcher.h"
#include "shapedb/molecules/RDMoleculeAnalytic.h"
#include "Orienter.h"
#include "shapedb/MemMapped.h"
#include "DatabaseStrutures.h"

using namespace std;

//searches a single directory
class FragmentSearcher
{
	GSSTreeSearcher gsstree; //for searching fragments, indexed by directory
	vector<DataIndex> indices; //keep these memory resident, gss trees index into here

	MemMapped sminaData;
	MemMapped molData;

public:
	//result of fragment search
	struct Result
	{
		unsigned pos;
		double score;

		Result(): pos(0), score(0) {}

		Result(unsigned p, double s): pos(p), score(s) {}
	};

	FragmentSearcher() {}
	 ~FragmentSearcher() {}

	//load in an existing index, return false on error
	bool read(const boost::filesystem::path& indirs);

	unsigned size() const { return indices.size(); }

	void search(GSSTreeSearcher::ObjectTree small, GSSTreeSearcher::ObjectTree big, vector< Result>& results); //TODO pharma

	//write sdf at pos to out
	void writeSDF(unsigned pos, ostream& out) const;

};

#endif /* FRAGMENTSEARCHER_H_ */
