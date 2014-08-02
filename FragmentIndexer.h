/*
 * FragmentIndexer.h
 *
 *  Created on: Jul 30, 2014
 *      Author: dkoes
 *
 *  This keeps track of molecular fragment data where the fragments are
 *  pieces of a reaction.  Fragments are not necessarily connected.
 *  It avoids storing redundant fragments (according to some rmsd cutoff).
 *  Fragments are assumed to be properly aligned to the corresponding scaffold.
 */

#ifndef FRAGMENTINDEXER_H_
#define FRAGMENTINDEXER_H_

#include <GraphMol/ROMol.h>
#include <GraphMol/Conformer.h>
#include <boost/unordered_map.hpp>
#include <vector>
#include <string>

#include "Orienter.h"
#include "MolMatcher.h"

using namespace std;

class FragmentIndexer
{
	struct Fragment
	{
		MolMatcher matcher;
		RDKit::ROMOL_SPTR frag; //connection table, etc
		vector<ECoords> coordinates; //store coordinates separately

		void initialize(const RDKit::Conformer& conf);
		//add conformer if it is more than cutoff away from current set
		void add(const RDKit::Conformer& conf, double cutoffSq);
	};

	double rmsdCutoffSq; //avoid duplicating fragments that are more similar than this

	boost::unordered_map<string, unsigned> fragmentPos; //indexed by smiles
	vector< Fragment > fragments; //indexed by fragmentPos, we maintain one molecule with multiple conformers

public:
	FragmentIndexer(): rmsdCutoffSq(0) {}
	FragmentIndexer(double rc): rmsdCutoffSq(rc*rc) {}

	//load in an existing index
	void read(istream& in);
	//write out index into specified file
	void write(ostream& out);

	//adds fragment to index if isn't within rmsdCufoff of previously added fragment
	void add(const RDKit::Conformer& conf);

	unsigned numFragments() const { return fragments.size(); }
	unsigned numFragmentConformers() const;

};

#endif /* FRAGMENTINDEXER_H_ */
