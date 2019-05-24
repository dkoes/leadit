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
#include <unordered_map>
#include <boost/filesystem.hpp>
#include <vector>
#include <string>
#include "shapedb/GSSTreeCreator.h"
#include "shapedb/molecules/RDMoleculeAnalytic.h"
#include "Orienter.h"
#include "MolMatcher.h"
#include "DatabaseStrutures.h"
#include "Reaction.h"

using namespace std;

class FragmentIndexer
{
	struct Fragment
	{
		MolMatcher matcher;
		RDKit::ROMOL_SPTR frag; //connection table, etc
		vector<ECoords> coordinates; //store coordinates separately
		vector<FragBondInfo> fragbonds;
		void initialize(const RDKit::Conformer& conf, const vector<FragBondInfo>& fragbonds);
		//add conformer if it is more than cutoff away from current set
		void add(const RDKit::Conformer& conf, double cutoffSq);

		void read(istream& in);
		void write(ostream& out) const;
	};


	//for writing out molecular data in an easy to access format
	//this class is responsible for outputting one slice of the data
	//functions as an iterator for gss tree creation
	class Outputter
	{
		vector<Fragment> *fragments; //reference to fragments
		unsigned position; //position in fragments array
		unsigned confposition; //position within fragment conformers
		unsigned stride;
		boost::filesystem::path dir;
		RDMolecule current;
		bool valid;

		vector<DataIndex> indices;
		ofstream *sminaData;
		ofstream *molData;
		ofstream *rdmolData;
		ofstream *pharmaData;

		void readNext();
		void setCurrent();

		Outputter(const Outputter& rhs); //disable copy constructor since bad things will happen
	public:
		Outputter(): fragments(NULL), position(0), confposition(0), stride(1), valid(false), sminaData(NULL), molData(NULL), pharmaData(NULL) {}

		Outputter(vector<Fragment> *frags, const boost::filesystem::path& d, unsigned start, unsigned strd);

		~Outputter()
		{
			if(sminaData) delete sminaData;
			if(molData) delete molData;
			if(pharmaData) delete pharmaData;
			if(rdmolData) delete rdmolData;
			sminaData = molData = pharmaData = NULL;
		}

		operator bool() const
		{
			return valid;
		}

		//current mol
		const RDMolecule& operator*() const
		{
			return current;
		}

		void operator++()
		{
			readNext();
		}

		const boost::filesystem::path& getDir() const { return dir; }

		void finish();

	};


	double rmsdCutoffSq; //avoid duplicating fragments that are more similar than this

	typedef std::unordered_map<string, unsigned> FragmentMap;
	FragmentMap fragmentPos; //indexed by smiles
	vector< Fragment > fragments; //indexed by fragmentPos, we maintain one molecule with multiple conformers

public:
	FragmentIndexer(): rmsdCutoffSq(0) {}
	FragmentIndexer(double rc): rmsdCutoffSq(rc*rc) {}

	//load in an existing index, return false on error
	bool read(const vector<boost::filesystem::path>& indirs);
	//stripe out index into specified directories - this will overwrite everything with the current data
	void write(const vector<boost::filesystem::path>& outdirs);

	//adds fragment to index if isn't within rmsdCufoff of previously added fragment
	void add(const RDKit::Conformer& conf, const vector<FragBondInfo>& fragbonds);

	unsigned numFragments() const { return fragments.size(); }
	unsigned numFragmentConformers() const;

};

#endif /* FRAGMENTINDEXER_H_ */
