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
	typedef shared_ptr<ChemicalReaction> ChemRxnPtr;
	ChemRxnPtr reverseRxn; //reaction
	ROMOL_SPTR product; //this is actually the left of the reaction
	vector<ROMOL_SPTR> reactants;
	ROMOL_SPTR core;
	vector<int> productReactants; //what reactant each product atom belongs to
	boost::unordered_set<unsigned> coreScaffoldSet; //set of core indices in product
	vector<unsigned> connectingMapNums;

	void findCoreAndReactants();
	void uniqueCoresOnly(vector<MatchVectType>& matches);
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

	unsigned numPieces() const { return reactants.size(); }

	//use the reaction to break up the passed mol into its starting reactants and core scaffold
	//also compute the indices of the connecting atoms in each reactant and the core
	//return false on failure
	bool decompose(const ROMol& mol, vector<MOL_SPTR_VECT>& react,
			vector<ROMOL_SPTR>& core);

	friend ostream& operator<<(ostream &out, Reaction &r); //for debugging

	void write(ostream& out); //for serialization
	void read(istream& in);
};

#endif /* REACTION_H_ */
