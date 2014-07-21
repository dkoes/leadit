/*
 * ScaffoldIndexer.cpp
 *
 *  Created on: Jul 21, 2014
 *      Author: dkoes
 */

#include <ScaffoldIndexer.h>

using namespace boost;

//create a new, empty index
void ScaffoldIndexer::initialize(const Reaction& rxn, double rmsdCut, double connectCut)
{
	rmsdCutoff = rmsdCut;
	connectCutoff = connectCut;
	//figure out core indices

}

//load in an existing index
void ScaffoldIndexer::read(istream& in)
{

}

//write out index into specified directory
void ScaffoldIndexer::write(ostream& out)
{

}

//finds the best scaffold cluster for the passed coordinates of the core scaffold
//which are assumed to be normalized to a standard order
//the cluster index is put in idx
//if the best match does not meet the matching criteria, return false
bool ScaffoldIndexer::findBest(ROMOL_SPTR core, vector<unsigned>& idx) const
{

}

//add a new scaffold conformation as represented by coords, will
//only create a new cluster if necessary, returns the cluster index
unsigned ScaffoldIndexer::addScaffold(ROMOL_SPTR core,  const vector<unsigned>& connect)
{

}
