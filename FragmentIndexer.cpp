/*
 * FragmentIndexer.cpp
 *
 *  Created on: Jul 30, 2014
 *      Author: dkoes
 */

#include <FragmentIndexer.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/RDKitBase.h>
#include <boost/graph/isomorphism.hpp>

using namespace RDKit;
using namespace boost;

#define ATOM_MAP_NUM "molAtomMapNumber"

void FragmentIndexer::Fragment::initialize(RDKit::CONFORMER_SPTR conf)
{
	//copy into a mol without conformers
	ROMol& mol = conf->getOwningMol();

	frag = ROMOL_SPTR(new ROMol(mol, true)); //quick copy - no conformers
	//but need to copy properties "molAtomMapNumber"
	for (ROMol::AtomIterator itr = mol.beginAtoms(), end = mol.endAtoms();
			itr != end; ++itr)
	{
		Atom *a = *itr;
		if(a->hasProp(ATOM_MAP_NUM))
		{
			int mapnum = 0;
			a->getProp(ATOM_MAP_NUM, mapnum);
			frag->getAtomWithIdx(a->getIdx())->setProp(ATOM_MAP_NUM, mapnum);
		}
	}

	//convert coordinates
	const RDGeom::POINT3D_VECT& coords = conf->getPositions();
	assert(coords.size() == frag->getNumAtoms());
	ECoords newc(coords.size(),3);
	for(unsigned i = 0, n = coords.size(); i < n; i++)
	{
		newc.row(i) << coords[i].x, coords[i].y, coords[i].z;
	}
	coordinates.push_back(newc);
}

//produce a unique vertex id for ismorphism that includes the atomic num and degree
struct VertexID
{
	const ROMol& mol;

	VertexID(const ROMol& m): mol(m) {}

	unsigned operator()(MolGraph::vertex_descriptor i) const
	{
		const Atom *a = mol.getAtomWithIdx(i);
		unsigned anum = a->getAtomicNum();
		assert(anum < 256);
		unsigned degree = a->getDegree();
		return anum*256+degree;
	}
};
//add conformer if it is more than cutoff away from current set
void FragmentIndexer::Fragment::add(RDKit::CONFORMER_SPTR conf, double cutoff)
{
	//first come up with mapping between mol atoms and conf - this might be slow
	//and an area for further optimization
	ROMol& confmol = conf->getOwningMol();
	const MolGraph& confgraph = confmol.getTopology();
	const MolGraph& molgraph = frag->getTopology();

	assert(num_vertices(confgraph) == num_vertices(molgraph));
	vector<MolGraph::vertex_descriptor> mapping(num_vertices(confgraph), 0);
	VertexID molv(*frag);
	VertexID confv(confmol);

	property_map<MolGraph, vertex_index_t>::type imap = get(vertex_index,confgraph);
	isomorphism(confgraph, molgraph,
			make_iterator_property_map(mapping.begin(), imap, mapping[0])
			);
			//confv, molv, 256*num_vertices(confgraph));

}

void FragmentIndexer::add(RDKit::CONFORMER_SPTR conf)
{
	ROMol& mol = conf->getOwningMol();
	//smiles are canonical and so can be used for indexing
	string smile = MolToSmiles(mol);

	if(fragmentPos.count(smile))
	{
		//already exists
		unsigned pos = fragmentPos[smile];
		fragments[pos].add(conf, rmsdCutoff);
	}
	else //need to create new fragment
	{
		unsigned pos = fragments.size();
		fragments.push_back(Fragment());
		fragments.back().initialize(conf);
		fragmentPos[smile] = pos;
	}
}

