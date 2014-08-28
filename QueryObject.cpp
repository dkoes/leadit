/*
 * QueryObject.cpp
 *
 *  Created on: Aug 28, 2014
 *      Author: dkoes
 */

#include <QueryObject.h>

#include <GraphMol/Conformer.h>
#include <GraphMol/FileParsers/MolWriters.h>
#include "shapedb/molecules/RDMoleculeAnalytic.h"
#include "shapedb/MappableOctTree.h"
#include "DatabaseStrutures.h"

using namespace RDKit;
using namespace boost;

//apply orient to molecular coordinates of mol and then generate objecttree
//leaves coordinates untouched
GSSTreeSearcher::ObjectTree MolecularQueryObject::getObjectTree(
		const Orienter& orient) const
{
	MappableOctTree *objTree = NULL;
	if (mol.getNumAtoms() > 0)
	{ //nonempty
		const Conformer& conf = mol.getConformer();

		//compute transformed coordinates
		RDGeom::POINT3D_VECT coords = conf.getPositions();
		orient.reorient(coords);

		Conformer *newconf = new Conformer(coords.size());
		newconf->getPositions() = coords;

		ROMol *remol = new ROMol(mol, true); //quick copy does not copy conformations
		remol->addConformer(newconf); //takes ownership of memory

		//SDWriter debug("qref.sdf");
		//debug.write(*remol);

		RDMolecule rmol(DEFAULT_GSS_DIMENSION, DEFAULT_GSS_RESOLUTION,
		DEFAULT_GSS_PROBERADIUS);
		rmol.set(remol); //takes ownership

		objTree = MappableOctTree::create(
		DEFAULT_GSS_DIMENSION,
		DEFAULT_GSS_RESOLUTION, rmol);

		if (shrink > 0)
		{
			//reduce the object
			MGrid grid;
			objTree->makeGrid(grid, DEFAULT_GSS_RESOLUTION);
			grid.shrink(shrink);
			free(objTree);
			objTree = MappableOctTree::createFromGrid(grid);
		}
		else if (shrink < 0) //grow
		{
			//expand the object
			MGrid grid;
			objTree->makeGrid(grid, DEFAULT_GSS_RESOLUTION);
			grid.grow(-shrink);
			free(objTree);
			objTree = MappableOctTree::createFromGrid(grid);
		}
	}
	else //empty mol, create empty objtree
	{
		MGrid grid(DEFAULT_GSS_DIMENSION, DEFAULT_GSS_RESOLUTION);
		objTree = MappableOctTree::createFromGrid(grid);
	}

	if (invert) //treat as excluded vol
		objTree->invert();

	return shared_ptr<const MappableOctTree>(objTree, free);
}

