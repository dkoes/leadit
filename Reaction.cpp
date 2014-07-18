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
#include <lemon/list_graph.h>
#include <lemon/bellman_ford.h>
#include <rdkit/GraphMol/Subgraphs/SubgraphUtils.h>
#include <GraphMol/Substruct/SubstructMatch.h>

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
	findCore();
}

#define ATOM_MAP_NUM "molAtomMapNumber"
//identify the core of the reaction
void Reaction::findCore()
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

	//create a graph from the product
	using namespace lemon;
	ListGraph graph;
	unsigned N = product->getNumAtoms();
	graph.reserveNode(N);

	//do our own attribute handling, these are all indexed by vertex id
	ListGraph::NodeMap<int> react(graph, -1); //reactant atom belongs to, -1 if none
	ListGraph::NodeMap<bool> heavy(graph, false); //definitely heavy atom
	ListGraph::NodeMap<Atom*> atoms(graph, NULL);
	vector<ListGraph::Node> nodes;
	nodes.reserve(N);

	//create nodes
	for (ROMol::AtomIterator aitr = product->beginAtoms(),
			end = product->endAtoms(); aitr != end; ++aitr)
	{
		Atom *a = *aitr;
		ListGraph::Node n = graph.addNode();
		assert(nodes.size() == a->getIdx()); //must correctly map
		nodes.push_back(n);
		atoms[n] = a;
		heavy[n] = a->getAtomicNum() != 1; //this thinks [#1,#6] is H, but that's okay since H can't be part of the core
		if (a->hasProp(ATOM_MAP_NUM))
		{
			int mapnum = 0;
			a->getProp(ATOM_MAP_NUM, mapnum);
			if (reactantMapNumbers.count(mapnum))
			{
				react[n] = reactantMapNumbers[mapnum];
			}
		}
	}
	assert(N == nodes.size());
	//setup edges
	for (ROMol::BondIterator bitr = product->beginBonds(),
			end = product->endBonds(); bitr != end; ++bitr)
	{
		Bond *b = *bitr;
		graph.addEdge(nodes[b->getBeginAtomIdx()], nodes[b->getEndAtomIdx()]);
	}

	//identify the connecting heavy atoms for reactants
	vector<ListGraph::Node> connecting;
	vector<vector<Atom*> > connectAtoms(reactants.size()); //indexed by reactant pos

	for (unsigned i = 0; i < N; i++)
	{
		ListGraph::Node n = nodes[i];
		int r = react[n];
		if (r >= 0 && heavy[n])
		{
			for (ListGraph::IncEdgeIt e(graph, n); e != INVALID; ++e)
			{
				ListGraph::Node u = graph.u(e);
				if (u == n) //apparently doesn't standardize order
					u = graph.v(e);
				assert(u != n);
				if (react[u] != r && heavy[u])
				{
					//cross reactants with heavy atoms
					connecting.push_back(n);
					connectAtoms[r].push_back(atoms[u]); //core atoms
					break;
				}
			}
		}
	}

	coreScaffoldSet.clear();
	//compute the shortest paths between all these connecting atoms
	//the atoms along these paths make up the minimal scaffold
	typedef BellmanFord<ListGraph> BF;
	ListGraph::ArcMap<int> amap(graph, 1); //all edges equal
	for (unsigned i = 0, n = connecting.size(); i < n; i++)
	{
		ListGraph::Node s = connecting[i];
		coreScaffoldSet.insert(atoms[s]->getIdx());
		BF bf(graph, amap);
		bf.run(s);
		//for each path from s to member of connecting
		for (unsigned j = i + 1; j < n; j++)
		{
			ListGraph::Node t = connecting[j];
			while (t != s)
			{
				coreScaffoldSet.insert(atoms[t]->getIdx());
				t = bf.predNode(t);
			}
		}
	}

	//set the bonds between the scaffold atoms
	std::map<int, int> mapping;
	coreScaffold = vector<int>(coreScaffoldSet.begin(), coreScaffoldSet.end());
	core = ROMOL_SPTR(
			Subgraphs::pathToSubmol(*product,
					Subgraphs::bondListFromAtomList(*product, coreScaffold),
					false, mapping));

	//save the connecting atoms in core
	coreConnectAtoms.resize(connectAtoms.size()); //indexed first by reactant pos
	for (unsigned r = 0, nr = connectAtoms.size(); r < nr; r++)
	{
		for (unsigned i = 0, n = connectAtoms[r].size(); i < n; i++)
		{
			Atom *a = connectAtoms[r][i];
			if (mapping.count(a->getIdx()))
			{
				coreConnectAtoms[r].push_back(mapping[a->getIdx()]);
			}
		}
		assert(coreConnectAtoms[r].size() > 0);
	}
}

//break up mol, return false on failure
//there may be multiple possible results
bool Reaction::decompose(const ROMol& mol, vector<MOL_SPTR_VECT>& react,
		vector<ROMOL_SPTR>& core)
{
	ROMOL_SPTR m(MolOps::addHs(mol)); //for proper match must have hydrogens
	vector<ROMOL_SPTR> rvect; //the full molecule "product" is on the left
	rvect.push_back(m);
	react = reverseRxn->runReactants(rvect);
	core.clear();
	if (react.size() == 0)
		return false; //reaction failed

	//need to extract appropraite core for each result
	std::vector<MatchVectType> matchesHere;
	unsigned matchCount = 0;

	matchCount = SubstructMatch(*m, *product, matchesHere, false, true, false);

	if (matchCount != react.size())
	{
		//this probably will fail eventually, but I want to see the failure example
		//before I figure out the best way to fix it
		cerr << "Mismatch in run reactants/substruct match\n";
		abort();
	}
	for (unsigned i = 0, n = matchesHere.size(); i < n; i++)
	{
		const MatchVectType &match = matchesHere[i];
		vector<int> productAtoms;
		//each match is a pair (query, mol) where query is our product
		for (unsigned j = 0, num = match.size(); j < num; j++)
		{
			if (coreScaffoldSet.count(match[j].first))
			{
				productAtoms.push_back(match[j].second);
			}
		}
		ROMOL_SPTR thiscore = ROMOL_SPTR(
				Subgraphs::pathToSubmol(*m,
						Subgraphs::bondListFromAtomList(*m, productAtoms)));
		core.push_back(thiscore);
	}
	return true;
}

ostream& operator<<(ostream &out, Reaction &r)
{
	out << "Reaction: " << ChemicalReactionToRxnSmarts(*r.reverseRxn)
			<< "\n";
	out << "Product: " << MolToSmiles(*r.product) << "\n";
	out << "Core: " << MolToSmiles(*r.core) << "\n";
	return out;
}
