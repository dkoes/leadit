/*
 * DatabaseSearcher.cpp
 *
 *  Created on: Aug 26, 2014
 *      Author: dkoes
 */

#include <DatabaseSearcher.h>
#include <fstream>

DatabaseSearcher::DatabaseSearcher(const vector<filesystem::path>& dbs): dbpaths(dbs)
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

			fragments[f][p].resize(fragdirs.size());
			for(unsigned d = 0, nd = fragdirs.size(); d < nd; d++)
			{
				bool success = fragments[f][p][d].read(fragdirs[d]);
				if(!success)
					return;
			}
		}
	}

	valid = true;
}

unsigned long DatabaseSearcher::totalConformers() const
{
	unsigned long tot = 0;
	for(unsigned i = 0, n1 = fragments.size(); i < n1; i++)
	{ //for each scaffold
		unsigned long product = 1;
		for(unsigned j = 0, n2 = fragments[i].size(); j < n2; j++)
		{//for each reactant
			unsigned long sum = 0;
			for(unsigned k = 0, n3 = fragments[i][j].size(); k < n3; k++)
			{ //for each directory
				sum += fragments[i][j][k].size();
			}
			product *= sum;
		}
		tot += product;
	}
	return tot;
}



