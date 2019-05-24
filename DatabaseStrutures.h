/*
 * DatabaseStructures.h
 *
 *  Shared classes etc for Database routines.
 *  Created on: Aug 26, 2014
 *      Author: dkoes
 */

#ifndef DATABASESTRUCTURES_H_
#define DATABASESTRUCTURES_H_

#include <fstream>
#include <RDGeneral/StreamOps.h>
#include <GraphMol/RDKitBase.h>

//make these compile time constants for now
#define DEFAULT_GSS_DIMENSION (64)
#define DEFAULT_GSS_RESOLUTION (0.5)
#define DEFAULT_GSS_PROBERADIUS (1.4)

//in order to more extendable, search indices point to an
//index in a lookup table that provides the index for other data types
//this is what is located in the lookup table
struct DataIndex
{
	unsigned long molloc;
	unsigned long sminaloc;
	unsigned long pharmaloc;
	unsigned long rdmolloc;
	//pharma data?

	void write(std::ostream& out) const
	{
		RDKit::streamWrite(out, molloc);
		RDKit::streamWrite(out, sminaloc);
		RDKit::streamWrite(out, pharmaloc);
		RDKit::streamWrite(out, rdmolloc);
	}

	bool read(std::istream& in)
	{
		RDKit::streamRead(in, molloc);
		RDKit::streamRead(in, sminaloc);
		RDKit::streamRead(in, pharmaloc);
		RDKit::streamRead(in, rdmolloc);
		return (bool)in;
	}
};

struct FragBondInfo
{
	RDKit::Bond::BondType order; //bond order
	unsigned char idx; //in fragment
	unsigned char cmap; //connecting map num not in fragment

	FragBondInfo(): order(RDKit::Bond::UNSPECIFIED), idx(0), cmap(0) {}
	FragBondInfo(RDKit::Bond::BondType o, int i, int cm): order(o), idx(i), cmap(cm) {}
};

#endif /* DATAINDEX_H_ */
