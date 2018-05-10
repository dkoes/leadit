/*
 * FragmentSearcher.cpp
 *
 *  Created on: Aug 26, 2014
 *      Author: dkoes
 */

#include "FragmentSearcher.h"
#include "FragmentIndexer.h"
#include "shapedb/molecules/PMol.h"
#include <GraphMol/MolPickler.h>
#include <GraphMol/FileParsers/MolWriters.h>
#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/device/array.hpp>
using namespace boost;

bool FragmentSearcher::read(const boost::filesystem::path& dir)
{
	filesystem::path gss = dir / "gss";
	if(!gsstree.load(gss))
		return false;

	filesystem::path idx = dir / "indices";
	ifstream indxfile(idx.string().c_str());

	DataIndex di;
	while(di.read(indxfile))
	{
		indices.push_back(di);
	}

	//map in data files
	filesystem::path mdata = dir / "molData";
	if(!molData.map(mdata.string(), true, false))
		return false;

	filesystem::path rdata = dir / "rdmolData";
	if(!rdmolData.map(rdata.string(), true, false))
		return false;

	filesystem::path smdata = dir / "sminaData";
	sminaData.map(smdata.string(), true, false); //not checking return value because I haven't implemented this yet, hopefully I will remember to come back here and fix it after implementing it

	return true;
}

//perform a distance constraint search on our gss tree and return the results
//TODO: pharmacophore constraints as well
void FragmentSearcher::search(GSSTreeSearcher::ObjectTree small, GSSTreeSearcher::ObjectTree big, vector< Result>& results)
{
	results.clear();

	PositionResults res;
	cout << "smvol " << small->volume() << "\tbigvol " << big->volume() << "\n";
	gsstree.dc_search(small, big, GSSTreeSearcher::ObjectTree(), true, res);

	//okay, not the most efficient in the world..
	results.reserve(res.size());
	for(unsigned i = 0, n = res.size(); i < n; i++)
	{
		results.push_back(Result(res.getPosition(i), res.getScore(i)));
	}
}

ROMOL_SPTR FragmentSearcher::readRDMol(unsigned pos, vector<FragBondInfo>& fragbonds) const
{
	using namespace boost::iostreams;

	assert(pos < indices.size());
	const char *data = rdmolData.begin() + indices[pos].rdmolloc;
	const char *end = rdmolData.end();
	if(pos < indices.size()-1)
	{
		end = rdmolData.begin() + indices[pos+1].rdmolloc;
	}

	basic_array_source<char> input_source(data, end);
	stream<basic_array_source<char> > input_stream(input_source);

	ROMOL_SPTR frag = ROMOL_SPTR(new ROMol());
	MolPickler::molFromPickle(input_stream, *frag);

	unsigned char n = 0;
	streamRead(input_stream, n);
	fragbonds.resize(n);
	for(unsigned i = 0; i < n; i++)
	{
		streamRead(input_stream, fragbonds[i].order);
		streamRead(input_stream, fragbonds[i].idx);
		streamRead(input_stream, fragbonds[i].cmap);
	}

	return frag;
}

//return index of atom in m with mapnum, -1 otherwise
typedef int (*GetMapNum)(const Atom* a);
static int findAtomWithMapNum(const RDKit::ROMol& m,  int mapnum)
{
	for (ROMol::ConstAtomIterator itr = m.beginAtoms(), end = m.endAtoms(); itr != end; ++itr)
	{
		const Atom *mola = *itr;
		unsigned idx = mola->getIdx();
		int mnum = Reaction::getMapNum(mola);
		if(mnum == mapnum)
		{
			return idx;
		}
	}
	return -1;
}


void FragmentSearcher::writeSDF(ROMOL_SPTR fullmol, const Reaction::Decomposition& decomp, unsigned whichr, unsigned pos, const Orienter& orient, ostream& out) const
{
	vector<FragBondInfo> fragbonds;
	ROMOL_SPTR frag = readRDMol(pos, fragbonds);
	//reorient fragment
	Conformer& conf = frag->getConformer(); //assume single conf

	orient.unorient(conf.getPositions());

	//reconstruct full molecule
	RWMol mol(*decomp.removePiece(fullmol, whichr));
	//not a fan of the loos coupling, but new indices of frag should be oldi+orign after insertion
	int orign = mol.getNumAtoms();
	mol.insertMol(*frag);

	//create bonds
	for(unsigned i = 0, n = fragbonds.size(); i < n; i++)
	{
		//TODO - make this more efficient, shouldn't have to hunt for the right atom
		const FragBondInfo& f = fragbonds[i];
		int srci = findAtomWithMapNum(mol, f.cmap); assert(srci >= 0);
		int dsti = orign + f.idx;
		assert(mol.getAtomWithIdx(dsti) != NULL);

		if (mol.getBondBetweenAtoms(srci, dsti) == NULL && srci != dsti)
		{
			mol.addBond(srci, dsti, f.order);
		}
	}

	//write it out
	RDKit::SDWriter writer(&out);
	writer.setKekulize(false);
	writer.write(mol);
	writer.close(); //flushes
}
