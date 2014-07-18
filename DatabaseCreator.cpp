/*
 * DatabaseCreator.cpp
 *
 *  Created on: Jul 16, 2014
 *      Author: dkoes
 */

#include <DatabaseCreator.h>

//creation
DatabaseCreator::DatabaseCreator(const filesystem::path& dbpath,
		const Reaction& react, const DatabaseConfiguration& config) :
		rxn(react)
{

}

//addition
DatabaseCreator::DatabaseCreator(const filesystem::path& dbpath)
{

}

DatabaseCreator::~DatabaseCreator()
{
}

//true if no problems setting up creation
bool DatabaseCreator::isValid() const
{
	return false;
}

//add conformers in molfile, only adds data, does not generate indices
void DatabaseCreator::add(const filesystem::path& molfile)
{

}

//generate indices and write all data to disk
void DatabaseCreator::finalize()
{

}
