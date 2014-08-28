/*
 * pharmacophores.cpp
 *
 *  Created on: Aug 5, 2014
 *      Author: dkoes
 *
 *  Implementation of pharmacophore rdkit matching.
 */

#include "pharmacophores.h"
#include <boost/assign.hpp>
#include <rdkit/GraphMol/SmilesParse/SmilesParse.h>
#include <rdkit/GraphMol/Substruct/SubstructMatch.h>
#include <rdkit/RDGeneral/StreamOps.h>

using namespace RDKit;
using namespace boost;
using namespace std;

// pharmacophore definitions - taken from pharmer
const char *aromatic[] =
		{ "a1aaaaa1", "a1aaaa1", NULL };

const char * hydrogen_donor[] =
		{ "[#7!H0&!$(N-[SX4](=O)(=O)[CX4](F)(F)F)]", "[#8!H0&!$([OH][C,S,P]=O)]",
				"[#16!H0]", NULL };

const char * hydrogen_acceptor[] =
		{ "[#7&!$([nX3])&!$([NX3]-*=[!#6])&!$([NX3]-[a])&!$([NX4])&!$(N=C([C,N])N)]",
				"[$([O])&!$([OX2](C)C=O)&!$(*(~a)~a)]", NULL };

const char * positive_ion[] =
		{ "[+,+2,+3,+4]",
				//amidine
				"[$(CC)](=N)N",
				//guanidine
				"[$(C(N)(N)=N)]", "[$(n1cc[nH]c1)]", NULL };

const char * negative_ion[] =
		{ "[-,-2,-3,-4]", "C(=O)[O-,OH,OX1]", "[$([S,P](=O)[O-,OH,OX1])]",
				"c1[nH1]nnn1", "c1nn[nH1]n1", "C(=O)N[OH1,O-,OX1]", "C(=O)N[OH1,O-]",
				"CO(=N[OH1,O-])",
				//trifluoromethyl sulfonamide
				"[$(N-[SX4](=O)(=O)[CX4](F)(F)F)]", NULL };

const char *hydrophobic[] =
		{
				"a1aaaaa1",
				"a1aaaa1",
				//branched terminals as one point
				"[$([CH3X4,CH2X3,CH1X2,F,Cl,Br,I])&!$(**[CH3X4,CH2X3,CH1X2,F,Cl,Br,I])]",
				"[$(*([CH3X4,CH2X3,CH1X2,F,Cl,Br,I])[CH3X4,CH2X3,CH1X2,F,Cl,Br,I])&!$(*([CH3X4,CH2X3,CH1X2,F,Cl,Br,I])([CH3X4,CH2X3,CH1X2,F,Cl,Br,I])[CH3X4,CH2X3,CH1X2,F,Cl,Br,I])]([CH3X4,CH2X3,CH1X2,F,Cl,Br,I])[CH3X4,CH2X3,CH1X2,F,Cl,Br,I]",
				"*([CH3X4,CH2X3,CH1X2,F,Cl,Br,I])([CH3X4,CH2X3,CH1X2,F,Cl,Br,I])[CH3X4,CH2X3,CH1X2,F,Cl,Br,I]",
				//simple rings only; need to combine points to get good results for 3d structures
				"[C&r3]1~[C&r3]~[C&r3]1",
				"[C&r4]1~[C&r4]~[C&r4]~[C&r4]1",
				"[C&r5]1~[C&r5]~[C&r5]~[C&r5]~[C&r5]1",
				"[C&r6]1~[C&r6]~[C&r6]~[C&r6]~[C&r6]~[C&r6]1",
				"[C&r7]1~[C&r7]~[C&r7]~[C&r7]~[C&r7]~[C&r7]~[C&r7]1",
				"[C&r8]1~[C&r8]~[C&r8]~[C&r8]~[C&r8]~[C&r8]~[C&r8]~[C&r8]1",
				//aliphatic chains
				"[CH2X4,CH1X3,CH0X2]~[CH3X4,CH2X3,CH1X2,F,Cl,Br,I]",
				"[$([CH2X4,CH1X3,CH0X2]~[$([!#1]);!$([CH2X4,CH1X3,CH0X2])])]~[CH2X4,CH1X3,CH0X2]~[CH2X4,CH1X3,CH0X2]",
				"[$([CH2X4,CH1X3,CH0X2]~[CH2X4,CH1X3,CH0X2]~[$([CH2X4,CH1X3,CH0X2]~[$([!#1]);!$([CH2X4,CH1X3,CH0X2])])])]~[CH2X4,CH1X3,CH0X2]~[CH2X4,CH1X3,CH0X2]~[CH2X4,CH1X3,CH0X2]",
				// sulfur (apparently)
				"[$([S]~[#6])&!$(S~[!#6])]", NULL };

const std::vector<PharmacophoreFeature> pharmacophoreFeatures = assign::list_of
		(PharmacophoreFeature("PharmaTypeAromatic", aromatic))
		(PharmacophoreFeature("PharmaTypeHydrogenDonor", hydrogen_donor))
		(PharmacophoreFeature("PharmaTypeHydrogenAcceptor", hydrogen_acceptor))
		(PharmacophoreFeature("PharmaTypePositiveIon", positive_ion))
		(PharmacophoreFeature("PharmaTypeNegativeIon", negative_ion))
		(PharmacophoreFeature("PharmaTypeHydrophobic", hydrophobic));

#define NUMPHARMAS (6)

//initialize smarts
PharmacophoreFeature::PharmacophoreFeature(const string& n, const char **sm) :
		propName(n)
{
	if (sm != NULL)
	{
		while (*sm != NULL)
		{
			smarts.push_back(ROMOL_SPTR(SmartsToMol(*sm)));
			sm++;
		}
	}
}

void assignPharmacophoreAtomProperties(RDKit::ROMOL_SPTR mol)
{
	for (unsigned i = 0, n = pharmacophoreFeatures.size(); i < n; i++)
	{
		const PharmacophoreFeature& ph = pharmacophoreFeatures[i];
		for (unsigned s = 0, ns = ph.smarts.size(); s < ns; s++)
		{
			MatchVectType match;
			if (SubstructMatch(*mol, *ph.smarts[s], match))
			{
				//match some pharmacophore, label all the participating atoms
				for (unsigned a = 0, na = match.size(); a < na; a++)
				{
					unsigned idx = match[a].second; //mol index
					Atom *atom = mol->getAtomWithIdx(idx);
					if (!atom->hasProp(ph.propName))
					{
						atom->setProp(ph.propName, i);
					}
				}
			}
		}
	}
}

//return bitmask of pharmacophore annotated on atom
unsigned atomPharmacophoreProps(RDKit::Atom *atom)
{
	unsigned ret = 0;
	for (unsigned i = 0, n = pharmacophoreFeatures.size(); i < n; i++)
	{
		const PharmacophoreFeature& ph = pharmacophoreFeatures[i];
		if (atom->hasProp(ph.propName))
			ret |= (1 << i);
	}
	return ret;
}

struct PharmaIndex
{
	struct Point
	{
		float x, y, z;
		Point() :
				x(0), y(0), z(0)
		{
		}

		Point(float X, float Y, float Z): x(X), y(Y), z(Z) {}
	};
	unsigned char start[NUMPHARMAS + 1];
	Point points[];
};

void writePharmacophoreProps(
		const vector<pair<RDGeom::Point3D, unsigned> >& props, ostream& out)
{
	//sort by pharma type
	vector< vector<PharmaIndex::Point> > pharmas;
	pharmas.resize(NUMPHARMAS);

	for(unsigned i = 0, n = props.size(); i < n; i++)
	{
		const RDGeom::Point3D& pt = props[i].first;
		unsigned mask = props[i].second;
		unsigned phpos = ffs(mask)-1; //least sig position is 1
		assert(phpos < NUMPHARMAS);
		pharmas[phpos].push_back(PharmaIndex::Point(pt.x, pt.y, pt.z));
	}

	//output the starting point of each pharma type
	unsigned pos = 0;
	for(unsigned i = 0, n = pharmas.size(); i < n; i++)
	{
		streamWrite(out, pos);
		pos += pharmas.size();
	}
	streamWrite(out, pos);
	//now output points
	for(unsigned i = 0, n = pharmas.size(); i < n; i++)
	{
		for(unsigned j = 0, m = pharmas[i].size(); j < m; j++)
		{
			const PharmaIndex::Point& pt = pharmas[i][j];
			streamWrite(out, pt.x);
			streamWrite(out, pt.y);
			streamWrite(out, pt.z);
		}
	}

}
