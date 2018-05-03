/*
 * FragmentSearcher.cpp
 *
 *  Created on: Aug 26, 2014
 *      Author: dkoes
 */

#include <FragmentSearcher.h>
#include "shapedb/molecules/PMol.h"

using namespace boost;

bool FragmentSearcher::read(const boost::filesystem::path& dir)
{
	filesystem::path gss = dir / "gss";
	if(!gsstree.load(gss))
		return false;

	filesystem::path idx = dir / "indices";
	ifstream indxfile(idx.string().c_str());

	DataIndex di;
	while(di.read(indxfile))
	{
		indices.push_back(di);
	}

	//map in data files
	filesystem::path mdata = dir / "molData";
	if(!molData.map(mdata.string(), true, false))
		return false;

	filesystem::path smdata = dir / "sminaData";
	sminaData.map(smdata.string(), true, false); //not checking return value because I haven't implemented this yet, hopefully I will remember to come back here and fix it after implementing it

	return true;
}

//perform a distance constraint search on our gss tree and return the results
//TODO: pharmacophore constraints as well
void FragmentSearcher::search(GSSTreeSearcher::ObjectTree small, GSSTreeSearcher::ObjectTree big, vector< Result>& results)
{
	results.clear();

	PositionResults res;
	cout << "smvol " << small->volume() << "\tbigvol " << big->volume() << "\n";
	gsstree.dc_search(small, big, GSSTreeSearcher::ObjectTree(), true, res);

	//okay, not the most efficient in the world..
	results.reserve(res.size());
	for(unsigned i = 0, n = res.size(); i < n; i++)
	{
		results.push_back(Result(res.getPosition(i), res.getScore(i)));
	}
}

void FragmentSearcher::writeSDF(unsigned pos, const Orienter& orient, ostream& out) const
{
	PMolReaderSimple reader;
	assert(pos < indices.size());
	const char *data = molData.begin() + indices[pos].molloc;
	PMol *pmol = reader.readPMol(data);
	pmol->writeSDF(out, orient);
}
