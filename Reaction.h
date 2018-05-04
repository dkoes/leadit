/*
 * Reaction.h
 *
 *  Created on: Jul 16, 2014
 *      Author: dkoes
 *
 *  Stores and manipulates reaction information.
 *  In addition to performing matching, can identify the core scaffold.
 */

#ifndef REACTION_H_
#define REACTION_H_

#include <vector>
#include <boost/filesystem.hpp>
#include <boost/unordered_set.hpp>
#include <RDGeneral/RDLog.h>
#include <RDGeneral/utils.h>
#include <GraphMol/RDKitBase.h>
#include <Geometry/point.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/ChemReactions/Reaction.h>
#include <GraphMol/ChemReactions/ReactionParser.h>

using namespace std;
using namespace boost;
using namespace RDKit;

#define ATOM_MAP_NUM "molAtomMapNumber"

class Reaction
{
public:
	//types
	typedef shared_ptr<ChemicalReaction> ChemRxnPtr;
	struct Connection
	{
		//record information on the junction between fragments
		int coreIndex;
		int reactantIndex; //indices are from the molecule provided to decompose
		int coreMap; // mappings from reaction
		int reactantMap;
		Bond::BondType order;

		Connection() :
				coreIndex(-1), reactantIndex(-1), coreMap(-1), reactantMap(-1), order(Bond::UNSPECIFIED)
		{
		}

		//so we can cannonicalize ordering of connections
		bool operator<(const Connection& rhs) const
				{
			//unsigned casts put negative numbers last
			if (coreMap == rhs.coreMap)
				return (unsigned) reactantMap < (unsigned) rhs.reactantMap;
			else
				return (unsigned) coreMap < (unsigned) rhs.coreMap;
		}

	};

	struct Decomposition
	{
		//store information for a single fragmentation
		vector<int> core; //atom indices of core
		vector<vector<int> > pieces; //of remaining fragments, indexed by reactant index
		vector<vector<Connection> > connections; //connections from core to fragments, indexed by reactant index, should be sorted
	};

private:
	ChemRxnPtr reverseRxn; //reaction
	ROMOL_SPTR product; //this is actually the left of the reaction
	vector<ROMOL_SPTR> reactants;
	ROMOL_SPTR core;
	vector<int> productReactants; //what reactant each product atom belongs to
	boost::unordered_set<unsigned> coreScaffoldSet; //set of core indices in product
	vector<unsigned> connectingMapNums;

	void findCoreAndReactants();
	void uniqueCoresOnly(vector<MatchVectType>& matches);
	bool decompFromMatch(ROMOL_SPTR m, const MatchVectType &match, Decomposition& decomp);

public:

	Reaction()
	{
	}
	Reaction(const filesystem::path& rxfile);
	Reaction(const Reaction& r) :
			reverseRxn(r.reverseRxn), product(r.product), reactants(
					r.reactants), core(r.core), productReactants(
					r.productReactants), coreScaffoldSet(
					r.coreScaffoldSet), connectingMapNums(r.connectingMapNums)
	{
	}
	~Reaction()
	{
	}

	bool isValid() const
	{
		return core && core->getNumAtoms() > 0;
	}

	unsigned coreSize() const
	{
		return core->getNumHeavyAtoms();
	}
	unsigned connectingSize() const
	{
		return connectingMapNums.size();
	}
	const vector<unsigned>& getConnectingMapNums() const
	{
		return connectingMapNums;
	}

	unsigned numPieces() const
	{
		return reactants.size();
	}

	//use the reaction to break up the passed mol into its starting reactants and core scaffold
	//also compute the indices of the connecting atoms in each reactant and the core
	//return false on failure
	bool decompose(ROMOL_SPTR m, vector<MOL_SPTR_VECT>& react,
			vector<ROMOL_SPTR>& core);

	//identify the decomposition of molecule m by this reaction using its
	//indices, without extracting submols.  Return true if successful
	//decomp length is number of matches
	bool decompose(ROMOL_SPTR m, vector<Decomposition>& decomp);

	//return sub-mol of m consisting of atoms specified in atomindices
	static ROMOL_SPTR extractMol(ROMOL_SPTR m, const vector<int>& atomindices);

	friend ostream& operator<<(ostream &out, const Reaction &r); //for debugging

	void write(ostream& out); //for serialization
	void read(istream& in);
};

#endif /* REACTION_H_ */
