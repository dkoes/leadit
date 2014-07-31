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

void FragmentIndexer::Fragment::initialize(const RDKit::Conformer& conf)
{
	//copy into a mol without conformers
	ROMol& mol = conf.getOwningMol();

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
	const RDGeom::POINT3D_VECT& coords = conf.getPositions();
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

	typedef unsigned result_type;
	typedef MolGraph::vertex_descriptor argument_type;

	VertexID(const ROMol& m): mol(m) {}

	result_type operator()(MolGraph::vertex_descriptor i) const
	{
		unsigned ret = 0;
		const Atom *a = mol.getAtomWithIdx(i);
		ret = a->getDegree();
		ret *= 256;
		ret += a->getAtomicNum();
		assert(a->getAtomicNum() < 256);
		ret *= 2;
		ret += a->getIsAromatic();
		return ret;
	}

	unsigned max() const {
		return mol.getNumAtoms()*512;
	}
};
//add conformer if it is more than cutoff away from current set
void FragmentIndexer::Fragment::add(const Conformer& conf, double cutoffSq)
{
	//first come up with mapping between mol atoms and conf - this might be slow
	//and an area for further optimization
	ROMol& confmol = conf.getOwningMol();
	const MolGraph& confgraph = confmol.getTopology();
	const MolGraph& molgraph = frag->getTopology();

	assert(num_vertices(confgraph) == num_vertices(molgraph));
	vector<MolGraph::vertex_descriptor> mapping(num_vertices(confgraph), 0); //indexed by mol atomid to get conf atom id
	VertexID molv(*frag);
	VertexID confv(confmol);

	property_map<MolGraph, vertex_index_t>::type imap = get(vertex_index,confgraph);
	isomorphism(confgraph, molgraph,
			isomorphism_map(make_iterator_property_map(mapping.begin(), imap, mapping[0])).
			vertex_invariant1(confv).vertex_invariant2(molv));

	//copy coordinates from conf in the same order as mol
	ECoords confcoords(mapping.size(), 3);
	for(unsigned i = 0, n = mapping.size(); i < n; i++)
	{
		unsigned idx = mapping[i];
		const RDGeom::Point3D& pt = conf.getAtomPos(idx);
		confcoords.row(i) << pt.x,pt.y,pt.z;
	}

	//now see if there already exists a conformer that is within cutoff
	for(unsigned i = 0, n = coordinates.size(); i < n; i++)
	{
		double r = (coordinates[i]-confcoords).squaredNorm()/(double)confcoords.size();
		if(r < cutoffSq)
			return; //no need to add
	}

	coordinates.push_back(confcoords);
}

void FragmentIndexer::add(const Conformer& conf)
{
	ROMol& mol = conf.getOwningMol();
	//smiles are canonical and so can be used for indexing
	string smile = MolToSmiles(mol);

	if(fragmentPos.count(smile))
	{
		//already exists
		unsigned pos = fragmentPos[smile];
		fragments[pos].add(conf, rmsdCutoffSq);
	}
	else //need to create new fragment
	{
		unsigned pos = fragments.size();
		fragments.push_back(Fragment());
		fragments.back().initialize(conf);
		fragmentPos[smile] = pos;
	}
}

unsigned FragmentIndexer::numFragmentConformers() const
{
	unsigned total = 0;
	for(unsigned i = 0, n = fragments.size(); i < n; i++)
	{
		total += fragments[i].coordinates.size();
	}
	return total;
}

