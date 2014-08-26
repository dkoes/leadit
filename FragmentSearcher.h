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
#include "DataIndex.h"
using namespace std;

//searches a single directory
class FragmentSearcher
{
	GSSTreeSearcher gsstree; //for searching fragments, indexed by directory
	vector<DataIndex> indices; //keep these memory resident, gss trees index into here

	MemMapped sminaData;
	MemMapped molData;

public:
	FragmentSearcher() {}
	 ~FragmentSearcher() {}

	//load in an existing index, return false on error
	bool read(const boost::filesystem::path& indirs);

	unsigned size() const { return indices.size(); }
};

#endif /* FRAGMENTSEARCHER_H_ */
