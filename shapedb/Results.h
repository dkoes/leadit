/*
 * Results.h
 *
 *  Created on: Jun 3, 2013
 *      Author: dkoes
 *
 *  Base class for results containers. Default behavior is to do nothing at all.
 */

#ifndef RESULTS_H_
#define RESULTS_H_

#include <vector>
#include <string>

class Results
{
public:
	Results()
	{
	}
	virtual ~Results()
	{
	}

	virtual void clear()
	{
	}

	//add result that is located at position within the databegin array with score
	virtual void add(const char *databegin, unsigned position, double score) = 0;

	virtual void reserve(unsigned n)
	{
	}
	virtual unsigned size() const
	{
		return 0;
	}
};

//for ojects that just store a string identifier (null terminated)
class StringResults: public Results
{
	std::vector<std::string> strs;
	std::vector<double> scores;
	public:
	StringResults()
	{
	}
	virtual ~StringResults()
	{
	}

	virtual void clear()
	{
		strs.clear();
		scores.clear();
	}
	virtual void add(const char *data, unsigned position, double score)
	{
		const char * addr = data + position;
		strs.push_back(addr); //better be null terminated
		scores.push_back(score);
	}
	virtual void reserve(unsigned n)
	{
		strs.reserve(n);
		scores.reserve(n);
	}
	virtual unsigned size() const
	{
		return strs.size();
	}

	const string& getString(unsigned i) const { return strs[i]; }
	double getScore(unsigned i) const { return scores[i]; }
};

//for results where there aren't any actual objects, we just want the position
//stored in the tree
class PositionResults: public Results
{
	std::vector<unsigned> positions;
	std::vector<double> scores;
	public:
	PositionResults()
	{
	}
	virtual ~PositionResults()
	{
	}

	virtual void clear()
	{
		positions.clear();
		scores.clear();
	}
	virtual void add(const char *data, unsigned position, double score)
	{
		positions.push_back(position); //better be null terminated
		scores.push_back(score);
	}
	virtual void reserve(unsigned n)
	{
		positions.reserve(n);
		scores.reserve(n);
	}
	virtual unsigned size() const
	{
		return positions.size();
	}

	unsigned getPosition(unsigned i) const { return positions[i]; }
	double getScore(unsigned i) const { return scores[i]; }
};

#endif /* RESULTS_H_ */
