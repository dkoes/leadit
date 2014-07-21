/*
 * ScaffoldIndexer.h
 *
 *  Created on: Jul 21, 2014
 *      Author: dkoes
 *
 *  This class creates, maintains, and searches a scaffold index.
 *  The scaffold index keeps track of clusters of similar conformations
 *  of the core scaffold.  It is assumed to be small and all operations
 *  are in-memory.  If multiple dbpaths are provided, the index is mirrored
 *  across all slices for full redundancy.
 */

#ifndef SCAFFOLDINDEXER_H_
#define SCAFFOLDINDEXER_H_

#include <vector>
#include <boost/filesystem.hpp>
#include <rdkit/Geometry/point.h>
#include "Reaction.h"
using namespace std;

class ScaffoldIndexer
{
	double rmsdCutoff; //minimum minimized RMSD
	double connectCutoff; //minimum unminimized RMSD to connection points
	unsigned numAtoms; //number of atoms in core scaffold
	vector<unsigned> connectingAtoms; //indices of connectors

	struct ScaffoldInfo
	{
		vector<RDGeom::Point3D> center;
		unsigned count; //how many we've seen like this

		ScaffoldInfo(): count(0) {}

		void read(istream& in, unsigned N);
		void write(ostream& out);
	};

	vector<ScaffoldInfo> clusters;


public:

	ScaffoldIndexer(): rmsdCutoff(0), connectCutoff(0), numAtoms(0) {}
	~ScaffoldIndexer() {}

	//create a new, empty index
	void initialize(const Reaction& rxn, double rmsdCut, double connectCut);

	//load in an existing index
	void read(istream& in);
	//write out index into specified file
	void write(ostream& out);

	//finds the best scaffold cluster for the passed core scaffold
	//the cluster index is put in idx
	//if the best match does not meet the matching criteria, return false
	bool findBest(ROMOL_SPTR core, vector<unsigned>& idx) const;

	//add a new scaffold conformation as represented by coords, will
	//only create a new cluster if necessary, returns the cluster index
	unsigned addScaffold(ROMOL_SPTR core, const vector<unsigned>& connect);
};

#endif /* SCAFFOLDINDEXER_H_ */
