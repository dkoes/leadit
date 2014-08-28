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
	//pharma data?

	void write(std::ostream& out) const
	{
		RDKit::streamWrite(out, molloc);
		RDKit::streamWrite(out, sminaloc);
	}

	bool read(std::istream& in)
	{
		RDKit::streamRead(in, molloc);
		RDKit::streamRead(in, sminaloc);
		return in;
	}
};

#endif /* DATAINDEX_H_ */
