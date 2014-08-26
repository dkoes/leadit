/*
 * FragmentSearcher.cpp
 *
 *  Created on: Aug 26, 2014
 *      Author: dkoes
 */

#include <FragmentSearcher.h>

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
