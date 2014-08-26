/*
 * DataIndex.h
 *
 *  Created on: Aug 26, 2014
 *      Author: dkoes
 */

#ifndef DATAINDEX_H_
#define DATAINDEX_H_

#include <fstream>
#include <RDGeneral/StreamOps.h>

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
		streamWrite(out, molloc);
		streamWrite(out, sminaloc);
	}

	bool read(std::istream& in)
	{
		streamRead(in, molloc);
		streamRead(in, sminaloc);
		return in;
	}
};

#endif /* DATAINDEX_H_ */
