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
#include <GraphMol/ChemReactions/ReactionPickler.h>

using namespace std;
using namespace boost;
using namespace RDKit;

class Reaction
{
	typedef shared_ptr<ChemicalReaction> ChemRxnPtr;
	ChemRxnPtr reverseRxn; //reaction
	ROMOL_SPTR product;
	vector<ROMOL_SPTR> reactants;

	boost::unordered_set<unsigned>  coreScaffoldSet; //set of core indices in product
	vector<int> coreScaffold; //atom indices of core in product
	ROMOL_SPTR core;
	vector< vector<int> > coreConnectAtoms;

	void findCore();
public:
	Reaction() {}
	Reaction(const filesystem::path& rxfile);
	Reaction(const Reaction& r): reverseRxn(r.reverseRxn), product(r.product), reactants(r.reactants), core(r.core) {}
	~Reaction() {}

	bool isValid() const
	{
		return core && core->getNumAtoms() > 0;
	}

	//use the reaction to break up the passed mol into its starting reactants and core scaffold
	//return false on failure
	bool decompose(const ROMol& mol, vector<MOL_SPTR_VECT>& react, vector<ROMOL_SPTR>& core);
	friend ostream& operator<< (ostream &out, Reaction &r);
};

#endif /* REACTION_H_ */
