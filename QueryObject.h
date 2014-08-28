/*
 * QueryObject.h
 *
 *  Created on: Aug 28, 2014
 *      Author: dkoes
 */

#ifndef QUERYOBJECT_H_
#define QUERYOBJECT_H_

#include "shapedb/GSSTreeSearcher.h"
#include <GraphMol/RWMol.h>
#include "Orienter.h"

/* A query object holds some sort of volumizable data that can be reoriented.
 * Basically, I don't like the idea of apply rigid body transformations directly
 * to grid data, since I suspect this would exacerbate grid effects and it isn't
 * necessarily very efficient.  So instead, we apply the transformation to the
 * underlying data, and then produce the grid.
 *
 * Eventually I want to support grid-based data in addition to molecular data,
 * so generalize this with a base class
 */
class QueryObject
{
public:
	QueryObject() {}
	virtual ~QueryObject() {}

	//return search object after applying passed transformation
	virtual GSSTreeSearcher::ObjectTree getObjectTree(const Orienter& orient) const = 0;
};


class MolecularQueryObject: public QueryObject
{
	RDKit::RWMol mol;
	double shrink; //amount to reduce mol surface by
	bool invert; //true if should compute negative shape (e.g. receptor)
public:
	MolecularQueryObject(): shrink(0), invert(false) {}
	MolecularQueryObject(bool i): invert(i) {}
	MolecularQueryObject(const MolecularQueryObject& rhs): mol(rhs.mol), shrink(rhs.shrink), invert(rhs.invert) {}
	MolecularQueryObject(const RDKit::ROMol& m, double s, bool i): mol(m), shrink(s), invert(i) {}
	~MolecularQueryObject() {}

	MolecularQueryObject& operator=(const MolecularQueryObject& rhs)
	{
		mol = rhs.mol;
		shrink = rhs.shrink;
		invert = rhs.invert;
		return *this;
	}

	GSSTreeSearcher::ObjectTree getObjectTree(const Orienter& orient) const;

};
#endif /* QUERYOBJECT_H_ */
