/*
 * DatabaseSearcher.cpp
 *
 *  Created on: Aug 26, 2014
 *      Author: dkoes
 */

#include <DatabaseSearcher.h>
#include <fstream>
#include <rdkit/GraphMol/FileParsers/MolWriters.h>
#include "CommandLine2/CommandLine.h"

extern cl::opt<bool> Verbose;

DatabaseSearcher::DatabaseSearcher(const vector<filesystem::path>& dbs) :
		dbpaths(dbs)
{
	assert(dbpaths.size() > 0);

	filesystem::path rxnpath = dbpaths[0] / "rxninfo";
	ifstream rxnf(rxnpath.c_str());
	if (!rxnf)
	{
		cerr << "Could not read " << rxnpath << "\n";
		return;
	}
	rxn.read(rxnf);

	filesystem::path idxpath = dbpaths[0] / "scaffoldIdx";
	ifstream scaffin(idxpath.c_str());
	if (!scaffin)
	{
		cerr << "Could not read " << idxpath << "\n";
		return;
	}
	scaffoldIndex.read(scaffin);

	unsigned sz = scaffoldIndex.size();
	unsigned nr = rxn.numPieces();
	fragments.resize(sz);

	for (unsigned f = 0; f < sz; f++)
	{
		fragments[f].resize(nr);
		for (unsigned p = 0; p < nr; p++)
		{
			vector<filesystem::path> fragdirs;
			for (unsigned i = 0, n = dbpaths.size(); i < n; i++)
			{
				filesystem::path fragdir = dbpaths[i] / lexical_cast<string>(f)
						/ lexical_cast<string>(p);
				if (!filesystem::is_directory(fragdir))
				{
					cerr << "Could not read " << fragdir << "\n";
					return;
				}
				fragdirs.push_back(fragdir);
			}

			fragments[f][p].resize(fragdirs.size());
			for (unsigned d = 0, nd = fragdirs.size(); d < nd; d++)
			{
				bool success = fragments[f][p][d].read(fragdirs[d]);
				if (!success)
					return;
			}
		}
	}

	valid = true;
}

unsigned long DatabaseSearcher::totalConformers() const
{
	unsigned long tot = 0;
	for (unsigned i = 0, n1 = fragments.size(); i < n1; i++)
	{ //for each scaffold
		unsigned long product = 1;
		for (unsigned j = 0, n2 = fragments[i].size(); j < n2; j++)
		{ //for each reactant
			unsigned long sum = 0;
			for (unsigned k = 0, n3 = fragments[i][j].size(); k < n3; k++)
			{ //for each directory
				sum += fragments[i][j][k].size();
			}
			product *= sum;
		}
		tot += product;
	}
	return tot;
}

//given a reference compounds (which must match the reaction and be in the same coordinate
//system as the other query objects) and the desired reactant index to replace,
//identify the relevent scaffolds and search for matches that fit between small and big
void DatabaseSearcher::search(ROMOL_SPTR ref, unsigned reactant,
		const MolecularQueryObject& small, const MolecularQueryObject& origbig,
		Results& results)
{
	results.clear();
	//analyze reference compound, does it match?
	ROMOL_SPTR m(MolOps::addHs(*ref,false,true)); //for proper match must have hydrogens
	results.refmol = m;

	rxn.decompose(m, results.decomps);
	const Conformer& conf = m->getConformer();

	if(Verbose)
	{
		cout << results.decomps.size() << " decompositions\n";
	}
	//for each core conformer
	for (unsigned c = 0, nc = results.decomps.size(); c < nc; c++)
	{
		Orienter coreorient;
		const Reaction::Decomposition& decomp = results.decomps[c];

		//add rest of molecule to big (excluded region)
		MolecularQueryObject big(origbig, *decomp.removePiece(m, reactant, false)); //don't include core
		//don't want to include rest of molecule (including core) in inclusive shape
		MolecularQueryObject smallex(*decomp.removePiece(m, reactant), small.getShrink(), small.getInvert());

		//canonicalize coordinates and add translation to orientation
		ECoords coords;
		scaffoldIndex.createCanonicalCoords(conf, decomp, coords, coreorient);

		//find any matching scaffolds
		vector<unsigned> scaffolds;
		scaffoldIndex.findBest(coords, scaffolds); //should we check return value, or just go with best anyway?

		if(Verbose)
		{
			cout << "\t" << scaffolds.size() << " scaffolds\n";
		}
		//for each match
		BOOST_FOREACH(unsigned& s, scaffolds)
		{
			//compute necessary rotation for this scaffold
			Orienter scaffoldorient(coreorient);
			scaffoldIndex.addRotation(s, coords, scaffoldorient);
			results.orienters.push_back(scaffoldorient);

			//apply orientation to query structures
			//remove rest of mol from small
			GSSTreeSearcher::ObjectTree smallextree = smallex.getObjectTree(scaffoldorient);
			smallextree->invert();
			GSSTreeSearcher::ObjectTree inctree = small.getObjectTree(scaffoldorient);
			GSSTreeSearcher::ObjectTree miv = QueryObject::getIntersection(inctree, smallextree);
			GSSTreeSearcher::ObjectTree msv = big.getObjectTree(scaffoldorient);

			//do search
			//TODO: multi thread, although threads will need to be hoisted
			assert(reactant < fragments[s].size());
			vector<FragmentSearcher::Result> fragResults; //will store fragment indices
			for (unsigned d = 0, nd = fragments[s][reactant].size(); d < nd; d++) //over directories
			{
				FragmentSearcher& searcher = fragments[s][reactant][d];
				searcher.search(miv, msv, fragResults);
				//expand results to fully specify positions
				BOOST_FOREACH(FragmentSearcher::Result& fr, fragResults)
				{
					results.results.push_back(Result(s, reactant, d, results.orienters.size() - 1, c, fr));
				}
			}
		}
	}
}

void DatabaseSearcher::writeSDF(const Results& results, unsigned rindex, ostream& out) const
{
	const Result& res = results.results[rindex];
	const FragmentSearcher& searcher = fragments[res.scaffoldPos][res.reactantPos][res.dir];
	assert(res.orientIndex < results.orienters.size());
	const Orienter& orient = results.orienters[res.orientIndex];
	searcher.writeSDF(results.refmol, results.decomps[res.dindex], res.reactantPos, res.fragRes.pos, orient, out);
}


