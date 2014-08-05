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
#include <boost/filesystem.hpp>
#include <vector>
#include <string>
#include "shapedb/GSSTreeCreator.h"
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

	//in order to more extendable, search indices point to an
	//index in a lookup table that provides the index for other data types
	//this is what is located in the lookup table
	struct DataIndex
	{
		unsigned long rdloc;
		unsigned long coordloc;
		unsigned long sminaloc;
	};

	//for writing out molecular data in an easy to access format
	struct OutputDir
	{
		//the gss tree of the shapes


		//the rdkit molecular data

		//the coordinate data

		//smina format

	};
	double rmsdCutoffSq; //avoid duplicating fragments that are more similar than this

	boost::unordered_map<string, unsigned> fragmentPos; //indexed by smiles
	vector< Fragment > fragments; //indexed by fragmentPos, we maintain one molecule with multiple conformers

public:
	FragmentIndexer(): rmsdCutoffSq(0) {}
	FragmentIndexer(double rc): rmsdCutoffSq(rc*rc) {}

	//load in an existing index, return false on error
	bool read(const vector<boost::filesystem::path>& indirs);
	//stripe out index into specified directories - this will overwrite everything with the current data
	void write(const vector<boost::filesystem::path>& outdirs);

	//adds fragment to index if isn't within rmsdCufoff of previously added fragment
	void add(const RDKit::Conformer& conf);

	unsigned numFragments() const { return fragments.size(); }
	unsigned numFragmentConformers() const;

};

#endif /* FRAGMENTINDEXER_H_ */
