/*
 * MolMatcher.cpp
 *
 *  Created on: Jul 31, 2014
 *      Author: dkoes
 */

#include <MolMatcher.h>
#include <rdkit/GraphMol/SmilesParse/SmilesWrite.h>
using namespace RDKit;

void MolMatcher::initialize(RDKit::ROMOL_SPTR rmol, const vector<RDGeom::Point3D>& rcoords)
{
	refmol = rmol;
	refcoords = rcoords;
	//setup coords and matching arrays
	unsigned n = refmol->getNumAtoms();
	assert(n == refcoords.size());
	matching.clear();
	matching.resize(n, -1);
	refmatching.clear();
	refmatching.resize(n, -1);
	atomids.resize(n, 0);
	for (unsigned i = 0; i < n; i++)
		atomids[i] = i;
}


//return true if the passed atoms match
//recusively descends and sets matching, but only if true is returned
bool MolMatcher::isMatch(const RDKit::Conformer& testconf, const Atom *testatom, const Atom *refatom)
{
	unsigned tidx = testatom->getIdx();
	unsigned ridx = refatom->getIdx();
	if (matching[tidx] >= 0)
	{
		if (matching[tidx] == (int)ridx)
			return true; //already matched
		else
			return false; //already matched differently
	}
	if(refmatching[ridx] >= 0)
	{
		if(refmatching[ridx] == (int)tidx)
			return true; //should have been caught above
		else
			return false;
	}
	if (testatom->getAtomicNum() != refatom->getAtomicNum())
		return false;
	if (testatom->getIsAromatic() != refatom->getIsAromatic())
		return false;
	if (testatom->getDegree() != refatom->getDegree())
		return false;

	//need to save/restore everything since we can give up in the middle of neighbors
	vector<int> savedmatching = matching;
	vector<int> savedrefmatching = refmatching;
	//these atoms are compatible, must now consider the neighbors
	matching[tidx] = ridx;
	refmatching[ridx] = tidx;

	//set up possible matches in ref
	vector<unsigned> neigh_idx; neigh_idx.reserve(refatom->getDegree());
	const ROMol& refmol = refatom->getOwningMol();

	ROMol::ADJ_ITER nbrIdx, endNbrs;
	for (boost::tie(nbrIdx, endNbrs) = refmol.getAtomNeighbors(refatom);
			nbrIdx != endNbrs; ++nbrIdx)
	{
		unsigned nidx = *nbrIdx;
		neigh_idx.push_back(nidx);
	}

	const ROMol& testmol = testatom->getOwningMol();
	//now go thorugh the target atoms's neighors and try to find matches
	for (boost::tie(nbrIdx, endNbrs) = testmol.getAtomNeighbors(testatom);
			nbrIdx != endNbrs; ++nbrIdx)
	{
		unsigned nidx = *nbrIdx;
		int match = findMatch(testconf, testmol.getAtomWithIdx(nidx), neigh_idx);
		if(match < 0) //could not match
		{
			//undo this matching
			matching = savedmatching;
			refmatching = savedrefmatching;
			return false;
		}
	}

	return true; //everything matched
}

//find a match for testatom in testmol, this descends recursively
//it only looks for matches within ids
//return the index of the matching refmol atom, -1 if not matchable
//has the side effect of seting matching iff a match is made
int MolMatcher::findMatch(const RDKit::Conformer& testmol,
		const RDKit::Atom *testatom, const vector<unsigned>& ids)
{
	unsigned idx = testatom->getIdx();
	if (matching[idx] >= 0)
	{
		//have a match, but is it within aids?
		unsigned m = matching[idx];
		for(unsigned i = 0, n = ids.size(); i < n;i++)
		{
			if(ids[i] == m)
				return m;
		}
		return -1;
	}

	//have to find match

	for (unsigned i = 0, n = ids.size(); i < n; i++)
	{
		const Atom *refatom = refmol->getAtomWithIdx(ids[i]);
		if (isMatch(testmol, testatom, refatom))
		{
			return ids[i];
		}
	}
	//did not find any match, unset matching

	return -1;
}

//find a matching
bool MolMatcher::computeMatch(const RDKit::Conformer& test)
{
	//simplest check for matchability
	if (test.getNumAtoms() != refmol->getNumAtoms())
		return false;

	const ROMol& confmol = test.getOwningMol();

	//the common case (possibly only case) is that the indices directly map
	bool identical = true;
	for(unsigned i = 0, n = refmol->getNumAtoms(); i < n; i++)
	{
		const Atom *ratom = refmol->getAtomWithIdx(i);
		const Atom *catom = confmol.getAtomWithIdx(i);
		if(ratom->getAtomicNum() != catom->getAtomicNum() ||
				ratom->getDegree() != catom->getDegree() ||
				ratom->getIsAromatic() != catom->getIsAromatic())
		{
			identical = false;
			break;
		}
		matching[i] = i;
		refmatching[i] = i;
	}
	if(identical) //same atom types, now check bonds
	{
		ROMol::BondIterator ritr = refmol->beginBonds(); //no auto conversion to const iterator
		ROMol::BondIterator rend = refmol->endBonds();
		for (ROMol::ConstBondIterator citr = confmol.beginBonds(), cend = confmol.endBonds();
				ritr != rend && citr != cend; ++ritr, ++citr)
		{
			const Bond *rbond = *ritr;
			const Bond *cbond = *citr;
			if(rbond->getBeginAtomIdx() != cbond->getBeginAtomIdx() ||
					rbond->getEndAtomIdx() != cbond->getEndAtomIdx())
			{
				identical = false;
				break;
			}
		}
	}

	if(identical)
		return true; //don't have to do full graph matching

	/* 99.9% of the time, we won't get here, but occastionally the rdkit
	 * numers don't match up so we have to do a full match
	 */

	matching.assign(matching.size(), -1);
	refmatching.assign(refmatching.size(), -1);
	//must match all atoms (don't assume connected
	for (ROMol::ConstAtomIterator itr = confmol.beginAtoms(), end =
			confmol.endAtoms();
			itr != end; ++itr)
	{
		const Atom *a = *itr;
		int m = findMatch(test, a, atomids);
		if(m < 0)
			return false;
	}


	return true;
}
