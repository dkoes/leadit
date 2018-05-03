/*
 * Reaction.cpp
 *
 *  Created on: Jul 16, 2014
 *      Author: dkoes
 */

#include <Reaction.h>
#include <iostream>
#include <fstream>
#include <boost/unordered_map.hpp>
#include <boost/unordered_set.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include <rdkit/GraphMol/Subgraphs/SubgraphUtils.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/ChemReactions/ReactionPickler.h>
#include <GraphMol/MolPickler.h>
#include <RDGeneral/StreamOps.h>

Reaction::Reaction(const filesystem::path& rxfile)
{
	//parse out smarts reaction and identify the core
	ifstream smarts(rxfile.c_str());
	string smart;
	smarts >> smart;
	ChemicalReaction *r = RxnSmartsToChemicalReaction(smart);
	if (r == NULL)
		throw runtime_error("Could not parse reaction");
	reverseRxn = ChemRxnPtr(r);

	//be flexible about whether it is a forwards or backwards reaction
	//but canonicalize it to be product first
	if (reverseRxn->getNumProductTemplates() == 1) //forward
	{
		product = *reverseRxn->beginProductTemplates();
		reactants.reserve(reverseRxn->getNumReactantTemplates());
		for (MOL_SPTR_VECT::iterator itr = reverseRxn->beginReactantTemplates(),
				end = reverseRxn->endReactantTemplates(); itr != end; ++itr)
		{
			reactants.push_back(*itr);
		}
		//create reverse rxn
		ChemRxnPtr reverse = make_shared<ChemicalReaction>();
		reverse->addReactantTemplate(product);
		for (unsigned i = 0, n = reactants.size(); i < n; i++)
		{
			reverse->addProductTemplate(reactants[i]);
		}
		reverseRxn = reverse;
	}
	else if (reverseRxn->getNumReactantTemplates() == 1)
	{
		product = *reverseRxn->beginReactantTemplates();
		reactants.reserve(reverseRxn->getNumProductTemplates());
		for (MOL_SPTR_VECT::iterator itr = reverseRxn->beginProductTemplates(),
				end = reverseRxn->endProductTemplates(); itr != end; ++itr)
		{
			reactants.push_back(*itr);
		}
	}
	else
	{
		throw runtime_error("Invalid reaction - no single product.");
	}

	reverseRxn->initReactantMatchers();
	findCoreAndReactants();
}

//identify the core of the reaction and label the noncore with what reactant they belong to
void Reaction::findCoreAndReactants()
{
	//record what atom mapping number belongs to each reactant
	unordered_map<int, int> reactantMapNumbers;
	for (unsigned i = 0, n = reactants.size(); i < n; i++)
	{
		ROMOL_SPTR r = reactants[i];
		for (ROMol::AtomIterator itr = r->beginAtoms(), end = r->endAtoms();
				itr != end; ++itr)
		{
			Atom *a = *itr;
			if (a->hasProp(ATOM_MAP_NUM))
			{
				int mapnum = 0;
				a->getProp(ATOM_MAP_NUM, mapnum);
				reactantMapNumbers[mapnum] = i;
			}
		}
	}

	//do our own attribute handling, these are all indexed by atom idx
	const MolGraph& graph = product->getTopology();
	unsigned N = product->getNumAtoms();
	productReactants = vector<int>(N, -1); //reactant atom belongs to, -1 if none
	vector<bool> heavy(N, false); //definitely heavy atom

	//get attributes for nodes nodes
	for (ROMol::AtomIterator aitr = product->beginAtoms(),
			end = product->endAtoms(); aitr != end; ++aitr)
	{
		Atom *a = *aitr;
		unsigned idx = a->getIdx();
		heavy[idx] = a->getAtomicNum() != 1; //this thinks [#1,#6] is H, but that's okay since H can't be part of the core
		if (a->hasProp(ATOM_MAP_NUM))
		{
			int mapnum = 0;
			a->getProp(ATOM_MAP_NUM, mapnum);
			if (reactantMapNumbers.count(mapnum))
			{
				productReactants[idx] = reactantMapNumbers[mapnum];
			}
		}
	}

	//identify the connecting heavy atoms for reactants
	vector<vector<Atom*> > connectAtoms(reactants.size()); //indexed by reactant pos
	vector<unsigned> connecting; //the indices of the connecting atom in the product
	connectingMapNums.clear(); //the map numbers of connecting atoms

	for (ROMol::AtomIterator aitr = product->beginAtoms(),
			end = product->endAtoms(); aitr != end; ++aitr)
	{
		Atom *a = *aitr;
		unsigned idx = a->getIdx();
		int r = productReactants[idx];
		if (r >= 0 && heavy[idx])
		{
			ROMol::ADJ_ITER nbrIdx, endNbrs;
			for (boost::tie(nbrIdx, endNbrs) = product->getAtomNeighbors(a);
					nbrIdx != endNbrs; ++nbrIdx)
			{
				unsigned u = *nbrIdx;
				if (productReactants[u] != r && heavy[u])
				{
					//cross reactants with heavy atoms
					connecting.push_back(idx);
					assert(a->hasProp(ATOM_MAP_NUM));
					int mapnum = 0;
					a->getProp(ATOM_MAP_NUM, mapnum);
					connectingMapNums.push_back(mapnum);
					connectAtoms[r].push_back(product->getAtomWithIdx(u)); //core atoms
					break;
				}
			}

		}
	}

	coreScaffoldSet.clear();
	//compute the shortest paths between all these connecting atoms
	//the atoms along these paths make up the minimal scaffold
	assert(num_vertices(graph) == N);
	for (unsigned i = 0, n = connecting.size(); i < n; i++)
	{
		unsigned s = connecting[i];
		coreScaffoldSet.insert(s);
		vector<MolGraph::vertex_descriptor> p(N);
		MolGraph::vertex_descriptor source(s);
		p[source] = source;
		breadth_first_search(graph, source, visitor(
				make_bfs_visitor(
						record_predecessors(&p[0], on_tree_edge()))));

		//for each path from s to member of connecting
		for (unsigned j = i + 1; j < n; j++)
		{
			unsigned t = connecting[j];
			MolGraph::vertex_descriptor node(t);
			while (node != source)
			{
				coreScaffoldSet.insert(node);
				node = p[node];
			}
		}
	}

	//set the bonds between the scaffold atoms
	vector<int> coreScaffold(coreScaffoldSet.begin(), coreScaffoldSet.end());
	core = ROMOL_SPTR(
			Subgraphs::atomsToSubmol(*product, coreScaffold));
}

//figure out which reaction atom idx belongs to
//this isn't very efficient - does a depth first traversal
//but it does update molreactants to avoid some recomputation
static int getMolReact(ROMOL_SPTR m, unsigned idx, vector<int>& molreactants, vector<bool>& molincore)
{
	if(molreactants[idx] >= 0)
		return molreactants[idx]; //basecase
	if(molreactants[idx] == -3) //currently being visited, avoid loop
		return -1;
	if(molincore[idx])
		return -1;
	int old = molreactants[idx];
	molreactants[idx] = -3;
	//otherwise we are whatever our neighbor is
	ROMol::ADJ_ITER nbrIdx, endNbrs;
	Atom *a = m->getAtomWithIdx(idx);
	for (boost::tie(nbrIdx, endNbrs) = m->getAtomNeighbors(a);
			nbrIdx != endNbrs; ++nbrIdx)
	{
		int neighReact = getMolReact(m, *nbrIdx, molreactants, molincore);
		//TODO: should check that reactants aren't overlapping
		if(neighReact >= 0)
		{
			molreactants[idx] = neighReact;
			return neighReact;
		}
	}
	//made it all the way through without finding something
	molreactants[idx] = old; //should be -2
	return old;
}

//copy all mapnums from src to dest
static void copyMapNums(ROMol& src, ROMol& dest, map<int,int>& mapping)
{
    for(INT_MAP_INT::const_iterator mapIt=mapping.begin();
          mapIt!=mapping.end();++mapIt) {
    	const Atom* srca = src.getAtomWithIdx(mapIt->first);
    	Atom *desta = dest.getAtomWithIdx(mapIt->second);

    	if(srca->hasProp(ATOM_MAP_NUM))
    	{
    		int mapnum = 0;
    		srca->getProp(ATOM_MAP_NUM, mapnum);
    		desta->setProp(ATOM_MAP_NUM, mapnum);
    	}
    }
}

//copy all mapnums from src to dest
static void copyMapNums(ROMol& src, ROMol& dest, const MatchVectType & mapping)
{
	for(unsigned i = 0, n = mapping.size(); i < n; i++)
	{
    	const Atom* srca = src.getAtomWithIdx(mapping[i].first);
    	Atom *desta = dest.getAtomWithIdx(mapping[i].second);

    	if(srca->hasProp(ATOM_MAP_NUM))
    	{
    		int mapnum = 0;
    		srca->getProp(ATOM_MAP_NUM, mapnum);
    		desta->setProp(ATOM_MAP_NUM, mapnum);
    	}
	}
}

//return mapnum or -1 if not set
static int getMapNum(const Atom* a)
{
	if(a->hasProp(ATOM_MAP_NUM))
	{
		int mapnum = 0;
		a->getProp(ATOM_MAP_NUM, mapnum);
		return mapnum;
	}
	return -1;
}

//reduce matches to only those that differ on the core scaffold
void Reaction::uniqueCoresOnly(vector<MatchVectType>& matches)
{
	vector<MatchVectType> newmatches; newmatches.reserve(matches.size());
	vector<MatchVectType> corematches; corematches.reserve(matches.size());

	for(unsigned m = 0, nm = matches.size(); m < nm; m++)
	{
		const MatchVectType &match = matches[m];

		MatchVectType corematch; corematch.reserve(match.size());
		for(unsigned i = 0, n = match.size(); i < n; i++)
		{
			unsigned prodatom = match[i].first;

			if (coreScaffoldSet.count(prodatom))
			{
				corematch.push_back(match[i]);
			}
		}

		if(find(corematches.begin(), corematches.end(), corematch) == corematches.end())
		{
			//avoid unnecessary copying with swap
			newmatches.push_back(MatchVectType());
			swap(newmatches.back(), matches[m]);
			corematches.push_back(MatchVectType());
			swap(corematches.back(),corematch);
		}
	}

	swap(matches,newmatches);
}

//use match to break m into pieces
bool Reaction::decompFromMatch(ROMOL_SPTR m, const MatchVectType &match,
		Decomposition& decomp)
{
	decomp.core.clear(); //molcoreAtoms
	decomp.pieces.clear();
	decomp.pieces.resize(reactants.size()); //reactAtoms

	vector<int> molreactants(m->getNumAtoms(), -2); //indexed by atom index, what reactant each atom is part of
	vector<bool> molincore(m->getNumAtoms(), false); //bitvector index by atom of core or not
	copyMapNums(*product, *m, match);

	//each match is a pair (query, mol) where query is our product
	for (unsigned j = 0, num = match.size(); j < num; j++)
	{
		unsigned prodatom = match[j].first;
		unsigned molatom = match[j].second;
		if (coreScaffoldSet.count(prodatom))
		{
			molincore[molatom] = true;
			decomp.core.push_back(molatom);
		}
		molreactants[molatom] = productReactants[prodatom];
	}

	//for each reactant, need to identify the part of m that came from that
	//reactant but that isn't part of the core scaffold
	for (ROMol::AtomIterator itr = m->beginAtoms(), end = m->endAtoms();
			itr != end; ++itr)
	{
		Atom *mola = *itr;
		unsigned idx = mola->getIdx();

		if (!molincore[idx])
		{
			int r = getMolReact(m, idx, molreactants, molincore);
			assert(r >= 0);
			assert(r < (int )decomp.pieces.size());
			if (mola->getAtomicNum() != 1) //ignore hydrogens
				decomp.pieces[r].push_back(idx);
		}
	}

	//each atom of mol is now labeled with what reactant it belongs to, although
	//some core atoms may still be unlabeled
	//need to identify the connections between non-core and core
	decomp.connections.resize(reactants.size()); //index by reactant first
	for (ROMol::BondIterator itr = m->beginBonds(), end = m->endBonds();
			itr != end; ++itr)
	{
		Bond *bond = *itr;
		unsigned a = bond->getBeginAtomIdx();
		unsigned b = bond->getEndAtomIdx();

		if (molincore[a] != molincore[b])
		{
			int core = -1, react = -1;
			if(molincore[a]) {
				core = a;
				react = b;
			}
			else {
				core = b;
				react = a;
			}

			Connection conn;
			conn.coreIndex = core;
			conn.reactantIndex = react;
			conn.coreMap = getMapNum(m->getAtomWithIdx(core));
			conn.reactantMap = getMapNum(m->getAtomWithIdx(react));
			conn.order = bond->getBondType();

			int r = molreactants[react];
			assert(r >= 0);
			decomp.connections[r].push_back(conn);

		}
		else if(!molincore[a] && !molincore[b])
		{
			if(molreactants[a] != molreactants[b])
			{
				//do not support reactants that connect with one another, must go through core
				return false;
			}
		}
	}

	return true;
}

bool Reaction::decompose(ROMOL_SPTR m, vector<Decomposition>& decomp)
{
	decomp.clear();

	//need to extract appropraite core for each result
	vector<MatchVectType> matchesHere;

	//we need to match the full product pattern
	SubstructMatch(*m, *product, matchesHere);

	uniqueCoresOnly(matchesHere);

	for (unsigned i = 0, n = matchesHere.size(); i < n; i++)
	{
		Decomposition d;
		const MatchVectType &match = matchesHere[i];
		if (decompFromMatch(m, match, d))
		{
			decomp.push_back(d);
		}
	}
	return decomp.size() > 0;
}

//break up mol, return false on failure
//there may be multiple possible results
//map nums are set in results
//assumes mol has hydrogens added
bool Reaction::decompose(ROMOL_SPTR m, vector<MOL_SPTR_VECT>& pieces,
		vector<ROMOL_SPTR>& core)
{
	core.clear();
	pieces.clear();

	//need to extract appropraite core for each result
	vector<MatchVectType> matchesHere;

	//we need to match the full product pattern
	SubstructMatch(*m, *product, matchesHere);

	uniqueCoresOnly(matchesHere);

	for (unsigned i = 0, n = matchesHere.size(); i < n; i++)
	{
		const MatchVectType &match = matchesHere[i];
		vector<int> molcoreAtoms;
		vector<int> molreactants(m->getNumAtoms(), -2);
		vector<bool> molincore(m->getNumAtoms(), false);
		copyMapNums(*product, *m, match);
		//each match is a pair (query, mol) where query is our product
		for (unsigned j = 0, num = match.size(); j < num; j++)
		{
			unsigned prodatom = match[j].first;
			unsigned molatom = match[j].second;
			if (coreScaffoldSet.count(prodatom))
			{
				molincore[molatom] = true;
				molcoreAtoms.push_back(molatom);
			}
			molreactants[molatom] = productReactants[prodatom];

		}

		//for each reactant, need to identify the part of m that came from that
		//reactant but that isn't part of the core scaffold
		vector< vector<int> > reactAtoms(reactants.size());
		for (ROMol::AtomIterator itr = m->beginAtoms(), end = m->endAtoms();
				itr != end; ++itr)
		{
			Atom *mola = *itr;
			unsigned idx = mola->getIdx();

			if(!molincore[idx])
			{
				int r = getMolReact(m, idx, molreactants, molincore);
				assert(r >= 0);
				assert(r < (int)reactAtoms.size());
				if(mola->getAtomicNum() != 1) //ignore hydrogens
					reactAtoms[r].push_back(idx);
			}

		}

		//each atom of mol is now labeled with what reactant it belongs to, although
		//some core atoms may still be unlabeled
		//need to identify the bonds between non-core and core
		vector< vector< pair<unsigned,unsigned> > > molcorebonds; //indexed first by reactant, pair goes from core to noncore
		molcorebonds.resize(reactants.size());
		for (ROMol::BondIterator itr = m->beginBonds(), end = m->endBonds();
				itr != end; ++itr)
		{
			Bond *bond = *itr;
			unsigned a = bond->getBeginAtomIdx();
			unsigned b = bond->getEndAtomIdx();

			if(molincore[a] && !molincore[b])
			{
				int r = molreactants[b];
				assert(r >= 0);
				molcorebonds[r].push_back(make_pair(a,b));

			}
			else if(!molincore[a] && molincore[b])
			{
				int r = molreactants[a];
				assert(r >= 0);
				molcorebonds[r].push_back(make_pair(b,a));
			}
		}


		//extract reactants
		MOL_SPTR_VECT reacts;
		for(unsigned r = 0, nr = reactAtoms.size(); r < nr; r++)
		{
			std::map<int, int> mapping;
			ROMOL_SPTR rreact = ROMOL_SPTR(
					Subgraphs::atomsToSubmol(*m, reactAtoms[r], mapping));

			copyMapNums(*m, *rreact, mapping);
			reacts.push_back(rreact);
		}

		pieces.push_back(reacts);

		//extract core
		std::map<int, int> mapping;
		ROMOL_SPTR thiscore = ROMOL_SPTR(
				Subgraphs::atomsToSubmol(*m, molcoreAtoms, mapping));
		copyMapNums(*m, *thiscore, mapping);

		//store indices of connecting atoms, this order should match that of the reactant order
		for(unsigned r = 0, nr = molcorebonds.size(); r < nr; r++)
		{
			for(unsigned b = 0, nb = molcorebonds[r].size(); b < nb; b++)
			{
				unsigned molcoreatom = molcorebonds[r][b].first;
				assert(mapping.count(molcoreatom));
			}
		}

		core.push_back(thiscore);
	}
	return true;
}

ostream& operator<<(ostream &out, const Reaction &r)
{
	out << "Reaction: " << ChemicalReactionToRxnSmarts(*r.reverseRxn)
			<< "\n";
	out << "Product: " << MolToSmiles(*r.product) << "\n";
	out << "Core: " << MolToSmiles(*r.core) << "\n";
	return out;
}

//serialize out
void Reaction::write(ostream& out)
{
	ReactionPickler::pickleReaction(*reverseRxn, out);

	MolPickler::pickleMol(*product,out);
	unsigned n = reactants.size();
	streamWrite(out, n);
	for(unsigned i = 0; i < n; i++)
	{
		MolPickler::pickleMol(*reactants[i],out);
	}
	MolPickler::pickleMol(*core, out);

	n = productReactants.size();
	streamWrite(out, n);
	for(unsigned i = 0; i < n; i++)
	{
		streamWrite(out,productReactants[i]);
	}
	n = coreScaffoldSet.size();
	streamWrite(out, n);
	for(unordered_set<unsigned>::iterator itr = coreScaffoldSet.begin(),
			end = coreScaffoldSet.end(); itr != end; ++itr)
	{
		unsigned val = *itr;
		streamWrite(out,val);
	}
}

//serialize in
void Reaction::read(istream& in)
{
	reverseRxn = make_shared<ChemicalReaction>();
	ReactionPickler::reactionFromPickle(in, *reverseRxn);

	product = ROMOL_SPTR(new ROMol());
	MolPickler::molFromPickle(in, *product);

	unsigned n = 0;
	streamRead(in, n);
	reactants.resize(n);

	for(unsigned i = 0; i < n; i++)
	{
		reactants[i] = ROMOL_SPTR(new ROMol());
		MolPickler::molFromPickle(in, *reactants[i]);
	}

	core = ROMOL_SPTR(new ROMol());
	MolPickler::molFromPickle(in, *core);

	streamRead(in, n);
	productReactants.resize(n);
	for(unsigned i = 0; i < n; i++)
	{
		streamRead(in, productReactants[i]);
	}

	streamRead(in, n);
	for(unsigned i = 0; i < n; i++)
	{
		unsigned val = 0;
		streamRead(in, val);
		coreScaffoldSet.insert(val);
	}
	assert(coreScaffoldSet.size() == n);
}
