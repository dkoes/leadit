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
#include <GraphMol/FileParsers/MolWriters.h>
#include <MolMatcher.h>
#include "shapedb/GSSTreeCreator.h"
#include "shapedb/packers/MatcherPacker.h"
#include "shapedb/KSamplePartitioner.h"
#include "shapedb/molecules/PMol.h"
#include "Orienter.h"
#include "pharmacophores.h"

using namespace RDKit;
using namespace boost;

#define ATOM_MAP_NUM "molAtomMapNumber"

void FragmentIndexer::Fragment::initialize(const RDKit::Conformer& conf)
{
	//copy into a mol without conformers
	ROMol& mol = conf.getOwningMol();

	frag = ROMOL_SPTR(new ROMol(mol, true)); //quick copy - no conformers

	matcher.initialize(frag, conf.getPositions());
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
	bool isiso = matcher.computeMatch(conf);
	if(!isiso)
	{
		const ROMol& refmol = matcher.getReference();
		const ROMol& confmol = conf.getOwningMol();
		cout << MolToSmiles(refmol) << "\n";
		cout << MolToSmiles(confmol) << "\n";

		unsigned n = refmol.getNumAtoms();
		for(unsigned i = 0; i < n; i++)
		{
			const Atom *a = refmol.getAtomWithIdx(i);
			cout << i << ":" << a->getSymbol() << "," << a->getDegree() << "," << a->getIsAromatic() << ",[";
			ROMol::ADJ_ITER nbrIdx, endNbrs;
			for (boost::tie(nbrIdx, endNbrs) = refmol.getAtomNeighbors(a);
					nbrIdx != endNbrs; ++nbrIdx)
			{
				cout << *nbrIdx << ",";
			}
			cout << "]  ";
		}
		cout << "\n";
		for(unsigned i = 0; i < n; i++)
		{
			const Atom *a = confmol.getAtomWithIdx(i);
			cout << i << ":" << a->getSymbol() << "," << a->getDegree() << "," << a->getIsAromatic() << ",[";
			ROMol::ADJ_ITER nbrIdx, endNbrs;
			for (boost::tie(nbrIdx, endNbrs) = confmol.getAtomNeighbors(a);
					nbrIdx != endNbrs; ++nbrIdx)
			{
				cout << *nbrIdx << ",";
			}
			cout << "]  ";
		}
		cout << "\n";
		isiso = matcher.computeMatch(conf);
	}
	assert(isiso);
	const vector<int>& matching = matcher.getMatching(); //maps from conf to refmol

	//copy coordinates from conf in the same order as mol
	ECoords confcoords(matching.size(), 3);
	for(unsigned i = 0, n = matching.size(); i < n; i++)
	{
		int refidx = matching[i];
		assert(refidx >= 0);
		const RDGeom::Point3D& pt = conf.getAtomPos(i);
		confcoords.row(refidx) << pt.x,pt.y,pt.z;
	}

	//now see if there already exists a conformer that is within cutoff
	for(unsigned i = 0, n = coordinates.size(); i < n; i++)
	{
		double r = (coordinates[i]-confcoords).squaredNorm()/(double)confcoords.size();
		if(r < cutoffSq)
		{
			return; //no need to add
		}
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

bool FragmentIndexer::read(const vector<boost::filesystem::path>& indirs)
{
	return false; //TODO - implement
}



//stripe out index into specified directories - this will overwrite everything with the current data
void FragmentIndexer::write(const vector<boost::filesystem::path>& outdirs)
{
	//create slices -- eventually multithread
	vector<Outputter> writers; writers.reserve(outdirs.size());
	for(unsigned i = 0, n = outdirs.size(); i < n; i++)
	{
		writers.push_back(Outputter(&fragments, outdirs[i], i, outdirs.size()));
	}

	//TODO: thread
	//these two are stateless
	MatcherPacker packer;
	KSamplePartitioner topdown;

	for(unsigned i = 0, n = writers.size(); i < n; i++)
	{
		GSSLevelCreator leveler(&topdown, &packer);
		GSSTreeCreator creator(&leveler);
		filesystem::path d = writers[i].getDir();

		if (!creator.create<RDMolecule, Outputter>(d / "gss", writers[i], writers[i].getDimension(), writers[i].getResolution(), true))
		{
			cerr << "Error creating gss database " << writers[i].getDir() << "\n";
			exit(-1);
		}

		//output indices
		writers[i].finish();
	}
}


FragmentIndexer::Outputter::Outputter(vector<Fragment> *frags, const boost::filesystem::path& d, unsigned start, unsigned strd):
	fragments(frags), position(start), confposition(0), stride(strd), dir(d), sminaData(NULL), molData(NULL)
{
	if(position < fragments->size() && confposition < (*fragments)[position].coordinates.size())
	{
		valid = true;
		filesystem::path sd = d / "sminaData";
		sminaData = new ofstream(sd.string().c_str());
		filesystem::path md = d / "molData";
		molData = new ofstream(md.string().c_str());

		setCurrent();
	}
	else
	{
		valid = false;
	}
}

//setup the next current fragment, outputting data as we go
void FragmentIndexer::Outputter::Outputter::readNext()
{
	confposition++; //next conformation
	if(confposition >= (*fragments)[position].coordinates.size())
	{
		//reset to next fragment
		confposition = 0;
		position += stride;
	}
	if(position < fragments->size())
	{
		//still valid
		setCurrent();
	}
	else
	{
		valid = false;
	}
}


void FragmentIndexer::Outputter::Outputter::setCurrent()
{
	//create a new romol with a single conformer
	const Fragment& frag = (*fragments)[position];
	ROMol *mol = new ROMol(*frag.frag);

	Conformer *conf = new Conformer(mol->getNumAtoms());
	//set atom positions
	ECoords coords = frag.coordinates[confposition];
	//also extract any pharmacphoric properties
	vector< pair<RDGeom::Point3D, unsigned> > props;
	for(unsigned i = 0, n = coords.rows(); i < n; i++)
	{
		RDGeom::Point3D pt;
		pt.x = coords(i,0);
		pt.y = coords(i,1);
		pt.z = coords(i,2);

		conf->setAtomPos(i, pt);

		Atom *atm = mol->getAtomWithIdx(i);
		unsigned pharmas = atomPharmacophoreProps(atm);
		if(pharmas != 0)
			props.push_back(make_pair(pt, pharmas));
	}
	mol->addConformer(conf); //takes ownershipt of conf

	current.set(mol); //takes ownershp of mol
	current.setProperties(props);

	//write out molecular data and save positions
	DataIndex idx;

	idx.sminaloc = sminaData->tellp();
	//TODO: generate and write sminaData

	idx.molloc = molData->tellp();
	PMolCreator pmol(*mol, true);
	pmol.writeBinary(*molData);

	indices.push_back(idx);
}

void FragmentIndexer::Outputter::Outputter::finish()
{
	//output indices
	filesystem::path indfile = dir / "indices";
	ofstream out(indfile.string().c_str());

	for(unsigned i = 0, n = indices.size(); i < n; i++)
	{
		out.write((char*)&indices[i],sizeof(indices[i]));
	}
}
