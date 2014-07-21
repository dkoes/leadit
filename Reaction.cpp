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

#define ATOM_MAP_NUM "molAtomMapNumber"
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
	std::map<int, int> mapping;
	vector<int> coreScaffold(coreScaffoldSet.begin(), coreScaffoldSet.end());
	core = ROMOL_SPTR(
			Subgraphs::pathToSubmol(*product,
					Subgraphs::bondListFromAtomList(*product, coreScaffold),
					false));

}

//enforce a canonical ordering of atoms
static bool canonicalAtomCompare(const Atom* l, const Atom *r)
{
	bool lprop = l->hasProp(ATOM_MAP_NUM);
	bool rprop = r->hasProp(ATOM_MAP_NUM);
	if(lprop && rprop)
	{
		int lnum = 0, rnum = 0;
		l->getProp(ATOM_MAP_NUM, lnum);
		r->getProp(ATOM_MAP_NUM, rnum);

		return lnum < rnum;
	}
	else if(lprop)
	{
		return true;
	}
	else if(rprop)
	{
		return false; //mapping go first
	}
	else
	{
		if(l->getAtomicNum() != r->getAtomicNum())
			return l->getAtomicNum() < r->getAtomicNum();
		return l->getIdx() < r->getIdx(); //hopefully this is sufficient
	}

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

//break up mol, return false on failure
//there may be multiple possible results
bool Reaction::decompose(const ROMol& mol, vector<MOL_SPTR_VECT>& react, vector< vector< vector<unsigned> > >& rconnect,
		vector<ROMOL_SPTR>& core, vector<vector<unsigned> >& coreConnect)
{
	ROMOL_SPTR m(MolOps::addHs(mol)); //for proper match must have hydrogens
	core.clear();
	coreConnect.clear();
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
		vector<int> molcoreAtoms;
		vector<int> molreactants(m->getNumAtoms(), -2);
		vector<bool> molincore(m->getNumAtoms(), false);
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


		rconnect.push_back(vector< vector<unsigned> >(reactAtoms.size()));
		//extract reactants
		MOL_SPTR_VECT reacts;
		for(unsigned r = 0, nr = reactAtoms.size(); r < nr; r++)
		{
			std::map<int, int> mapping;
			ROMOL_SPTR rreact = ROMOL_SPTR(
					Subgraphs::pathToSubmol(*m,
							Subgraphs::bondListFromAtomList(*m, reactAtoms[r]),
							false, mapping));


			reacts.push_back(rreact);

			//save the remaped connection points
			for(unsigned j = 0, nj = molcorebonds[r].size(); j < nj; j++)
			{
				unsigned molidx = molcorebonds[r][j].second;
				assert(mapping.count(molidx) > 0);
				unsigned newidx = mapping[molidx];
				rconnect.back()[r].push_back(newidx);
			}
		}

		//extract core
		std::map<int, int> mapping;
		ROMOL_SPTR thiscore = ROMOL_SPTR(
				Subgraphs::pathToSubmol(*m,
						Subgraphs::bondListFromAtomList(*m, molcoreAtoms),
						false, mapping));



		//store indices of connecting atoms, this order should match that of the reactant order
		vector<unsigned> corei;
		for(unsigned r = 0, nr = molcorebonds.size(); r < nr; r++)
		{
			for(unsigned b = 0, nb = molcorebonds[r].size(); b < nb; b++)
			{
				unsigned molcoreatom = molcorebonds[r][b].first;
				assert(mapping.count(molcoreatom));
				unsigned connectatom = mapping[molcoreatom];
				corei.push_back(connectatom);
			}
		}

		coreConnect.push_back(corei);
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

//serialize out
void Reaction::write(ostream& out)
{
	ReactionPickler::pickleReaction(*reverseRxn, out);
	MolPickler::pickleMol(*product,out);
	unsigned n = reactants.size();
	out.write((char*)&n, sizeof(unsigned));
	for(unsigned i = 0; i < n; i++)
	{
		MolPickler::pickleMol(*reactants[i],out);
	}
	MolPickler::pickleMol(*core, out);

	n = productReactants.size();
	out.write((char*)&n, sizeof(unsigned));
	for(unsigned i = 0; i < n; i++)
	{
		out.write((char*)&productReactants[i],sizeof(int));
	}
	n = coreScaffoldSet.size();
	out.write((char*)&n, sizeof(unsigned));
	for(unordered_set<unsigned>::iterator itr = coreScaffoldSet.begin(),
			end = coreScaffoldSet.end(); itr != end; ++itr)
	{
		unsigned val = *itr;
		out.write((char*)&val,sizeof(unsigned));
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
	in.read((char*)&n, sizeof(unsigned));
	reactants.resize(n);

	for(unsigned i = 0; i < n; i++)
	{
		reactants[i] = ROMOL_SPTR(new ROMol());
		MolPickler::molFromPickle(in, *reactants[i]);
	}

	core = ROMOL_SPTR(new ROMol());
	MolPickler::molFromPickle(in, *core);

	in.read((char*)&n, sizeof(unsigned));
	productReactants.resize(n);
	for(unsigned i = 0; i < n; i++)
	{
		in.read((char*)&productReactants[i],sizeof(int));
	}

	in.read((char*)&n, sizeof(unsigned));
	for(unsigned i = 0; i < n; i++)
	{
		unsigned val = 0;
		in.read((char*)&val, sizeof(unsigned));
		coreScaffoldSet.insert(val);
	}
	assert(coreScaffoldSet.size() == n);
}
