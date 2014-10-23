/*
 * MolMatcher.h
 *
 *  Created on: Jul 31, 2014
 *      Author: dkoes
 */

#ifndef MOLMATCHER_H_
#define MOLMATCHER_H_

#include <rdkit/GraphMol/RDKitBase.h>
#include <rdkit/GraphMol/ROMol.h>
#include <rdkit/Geometry/point.h>
#include <boost/unordered_map.hpp>
#include <vector>

using namespace std;

/* class to perform graph matching between two molecules.
 * Is initialized with the reference molecule.
 * Will figure out the atom correspondences and compute the rmsd between the
 * ref mol and a provided test mol.
 *
 * This is necessary because boost's isomorphism function can't deal with
 * disconnected graphs nor is it easy to break ties with positions
 *
 * This is designed around several assumptions that appropriate for leadit:
 * -the common case is that there is a one-to-one matching between indices
 * -the two fragments are aligned to close atoms probably match
 */
class MolMatcher
{
	RDKit::ROMOL_SPTR refmol;
	vector<RDGeom::Point3D> refcoords;
	vector<unsigned> atomids; //just the numbers from 0-n
	vector<int> matching; //map from test to ref
	vector<int> refmatching; //map from ref to test

	bool isMatch(const RDKit::Conformer& testconf, const RDKit::Atom *testatom, const RDKit::Atom *refatom);
	int findMatch(const RDKit::Conformer& testmol, const RDKit::Atom *testatom, const vector<unsigned>& ids);
public:

	MolMatcher() {}
	MolMatcher(RDKit::ROMOL_SPTR rmol, const vector<RDGeom::Point3D>& rcoords) { initialize(rmol, rcoords); }

	void initialize(RDKit::ROMOL_SPTR rmol, const vector<RDGeom::Point3D>& rcoords);

	//computes a correspondence between the ref mol and test (exhaustively)
	//and returns the rmsd; returns infinity if unmatchable
	bool computeMatch(const RDKit::Conformer& test);

	//provide access to match vector
	const vector<int>& getMatching() const { return matching; }

	const RDKit::ROMol& getReference() const { return *refmol; }

	void read(istream& in);
	void write(ostream& out) const;
};
#endif /* MOLMATCHER_H_ */
