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
#include <eigen3/Eigen/Core>
#include "Reaction.h"
#include "Orienter.h"
using namespace std;


class ScaffoldIndexer
{
	double rmsdCutoffSq; //minimum minimized RMSD squared
	double connectCutoffSq; //minimum unminimized RMSD to connection points squared
	unsigned numAtoms; //number of atoms in core scaffold
	unordered_set<unsigned> connectingMapNums; //map nums of connecting atoms

	struct ScaffoldInfo
	{
		ECoords center;
		unsigned count; //how many we've seen like this

		ScaffoldInfo(): count(0) {}

		void read(istream& in, unsigned N);
		void write(ostream& out);
	};

	vector<ScaffoldInfo> clusters;

	//calculate connecting and overall RMSD, return true if cutoffs are made
	bool calcRMSDSquares(const ECoords& a, const ECoords& b, double& connectRMSDSq, double& totalRMSDSq) const;

	//put coordinates of atom into canonical form, sorted by mapping id
	//and with connecting atoms first
	void createCanonicalCoords(const Conformer& core, ECoords& coords, Orienter& orient) const;

	static EMatrix3 computeRotation(const ECoords& ref, const ECoords& fit);
public:

	ScaffoldIndexer(): rmsdCutoffSq(0), connectCutoffSq(0), numAtoms(0) {}
	~ScaffoldIndexer() {}

	//create a new, empty index
	void initialize(const Reaction& rxn, double rmsdCut, double connectCut);

	//load in an existing index
	void read(istream& in);
	//write out index into specified file
	void write(ostream& out);

	//finds the best scaffold cluster for the passed coordinates
	//the cluster index is put in idx
	//if the best match does not meet the matching criteria, return false
	bool findBest(const ECoords& coords, vector<unsigned>& idx) const;

	//add a new scaffold conformation as represented by coords, will
	//only create a new cluster if necessary, returns the cluster index
	//returns the orienter needed to align molecule to chosen scaffold
	unsigned addScaffold(const Conformer& core, Orienter& orient);

	unsigned size() const { return clusters.size(); }

	//return number of molecules that have this scaffold conformation
	unsigned getCnt(unsigned i) const { return clusters[i].count; }

	void dumpCounts(ostream& out) const;
};

#endif /* SCAFFOLDINDEXER_H_ */
