/*
 * pharmacophores.h
 *
 *  Created on: Aug 5, 2014
 *      Author: dkoes
 *
 *  Routines for recognizing pharmacophore features in rdkit molecules
 */

#ifndef PHARMACOPHORES_H_
#define PHARMACOPHORES_H_

#include <rdkit/GraphMol/RDKitBase.h>
#include <rdkit/GraphMol/ROMol.h>
#include <vector>
#include <string>

struct PharmacophoreFeature
{
	std::string propName; //how stored in atoms
	std::vector<RDKit::ROMOL_SPTR> smarts;

	PharmacophoreFeature() {}
	PharmacophoreFeature(const std::string& n, const char **smarts);
};

extern const std::vector<PharmacophoreFeature> pharmacophoreFeatures; //stores names of pharmacophore properties

//assign pharmacophore types (possibly multiple) to each atom in mol
void assignPharmacophoreAtomProperties(RDKit::ROMOL_SPTR mol);

//return bitmask of any pharmacophore feature atom is annotated with, -1 if none
unsigned atomPharmacophoreProps(RDKit::Atom *atom);

#endif /* PHARMACOPHORES_H_ */
