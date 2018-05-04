/*
 * ScaffoldIndexer.cpp
 *
 *  Created on: Jul 21, 2014
 *      Author: dkoes
 */

#include <ScaffoldIndexer.h>
#include <eigen3/Eigen/Eigenvalues>
#include <eigen3/Eigen/SVD>
#include <RDGeneral/StreamOps.h>
using namespace Eigen;
using namespace boost;

//create a new, empty index
void ScaffoldIndexer::initialize(const Reaction& rxn, double rmsdCut,
		double connectCut)
{
	rmsdCutoffSq = rmsdCut * rmsdCut;
	connectCutoffSq = connectCut * connectCut;
	numAtoms = rxn.coreSize();
	connectingMapNums.insert(rxn.getConnectingMapNums().begin(),
			rxn.getConnectingMapNums().end());
}

//load in an existing index
void ScaffoldIndexer::read(istream& in)
{
	streamRead(in, rmsdCutoffSq);
	streamRead(in, connectCutoffSq);
	streamRead(in, numAtoms);

	unsigned n = 0;
	streamRead(in, n);
	for (unsigned i = 0; i < n; i++)
	{
		unsigned val = 0;
		streamRead(in, val);
		connectingMapNums.insert(val);
	}

	//number of clusters
	streamRead(in, n);
	clusters.resize(n);
	for(unsigned i = 0; i < n; i++)
	{
		clusters[i].read(in, numAtoms);
	}
}

//read in scaffold info, N is the number of atoms in a core
void ScaffoldIndexer::ScaffoldInfo::read(istream& in, unsigned int N)
{
	streamRead(in, count);
	center.resize(N, 3);
	for(unsigned i = 0, n = center.size(); i < n; i++)
	{
		float val = 0;
		streamRead(in, val);
		center(i) = val;
	}


}

void ScaffoldIndexer::ScaffoldInfo::write(ostream& out)
{
	streamWrite(out, count);
	for(unsigned i = 0, n = center.size(); i < n; i++)
	{
		float val = center(i);
		streamWrite(out, val);
	}
}

//write out index into specified directory
void ScaffoldIndexer::write(ostream& out)
{
	streamWrite(out, rmsdCutoffSq);
	streamWrite(out, connectCutoffSq);
	streamWrite(out, numAtoms);

	unsigned n = connectingMapNums.size();
	streamWrite(out, n);

	BOOST_FOREACH(unsigned val, connectingMapNums)
	{
		streamWrite(out, val);
	}

	n = clusters.size();
	streamWrite(out, n);
	for(unsigned i = 0; i < n; i++)
	{
		assert(clusters[i].center.rows() == numAtoms);
		clusters[i].write(out);
	}
}

static void dumpXYZ(const ECoords& c, const char *name)
{
	cout << c.rows() << "\n";
	cout << name << "\n";
	for (unsigned i = 0, nr = c.rows(); i < nr; i++)
	{
		cout << "C\t" << c.row(i) << "\n";
	}
}

//calculate rotation matrix to move b onto a, assumes they are already centered
EMatrix3 ScaffoldIndexer::computeRotation(const ECoords& a, const ECoords& b)
{
	//find the minimal all atom rmsd alignment using the Kabsch algorithm
	//both are assumed centered
	//compute the covariance matrix
	EMatrix3 A = a.transpose() * b;

	JacobiSVD<EMatrix3> svd(A, ComputeFullU | ComputeFullV);
	EMatrix3 U = svd.matrixU();
	EMatrix3 V = svd.matrixV();
	EMatrix3 R = U * V.transpose();
	double d = R.determinant();
	if (d < 0) //avoid reflections
	{
		V.col(2) *= -1;
		R = U * V.transpose(); //for a 3x3 matrix this flips the det
	}

	return R.transpose();
}

//calculate connecting and overall RMSD, return true if cutoffs are made
//coordinates are assumed to be canonical and therefore centered
bool ScaffoldIndexer::calcRMSDSquares(const ECoords& a, const ECoords& b,
		double& connectRMSDSq, double& totalRMSDSq) const
		{
	double connectSum = 0;
	double totalSum = 0;

	assert(a.rows() == b.rows());
	assert(a.rows() == numAtoms);
	unsigned nc = connectingMapNums.size();

	ECoords newb = b * computeRotation(a, b);

	//connecting atoms are defined to be first
	unsigned i;
	for (i = 0; i < nc; i++)
	{
		connectSum += (a.row(i) - newb.row(i)).squaredNorm();
	}
	totalSum = connectSum;
	for (/* i okay */; i < numAtoms; i++)
	{
		totalSum += (a.row(i) - newb.row(i)).squaredNorm();
	}

	connectRMSDSq = connectSum / nc;
	totalRMSDSq = totalSum / numAtoms;
	return (connectRMSDSq < connectCutoffSq) && (totalRMSDSq < rmsdCutoffSq);
}

//for sorting
struct MapNumInfo
{
	unsigned mapnum;
	unsigned idx;

	MapNumInfo() :
			mapnum(0), idx(0)
	{
	}
	MapNumInfo(unsigned i, unsigned m) :
			mapnum(m), idx(i)
	{
	}

	bool operator<(const MapNumInfo& rhs) const
			{
		return mapnum < rhs.mapnum;
	}

	bool operator==(const MapNumInfo& rhs) const
			{
		return mapnum == rhs.mapnum;
	}
};

//return true if most of the weight of the col is on the positive side
//for purposes of standardizing an alignment
static bool isPositiveBiased(const ECoords& coords, unsigned col)
{
	unsigned A = (col + 1) % 3;
	unsigned B = (col + 2) % 3;
	float total = 0;
	for (unsigned i = 0, nr = coords.rows(); i < nr; i++)
	{
		float val = coords(i, A) * coords(i, A) + coords(i, B) * coords(i, B);
		if (coords(i, col) < 0)
			total -= val;
		else
			total += val;
	}
	return total > 0;
}

//singleton for rotation matrices
struct RotationMatrices
{
	EMatrix3 xrot;
	EMatrix3 yrot;
	EMatrix3 zrot;

	RotationMatrices()
	{
		xrot << 1, 0, 0,
				0, -1, 0,
				0, 0, -1;
		yrot << -1, 0, 0,
				0, 1, 0,
				0, 0, -1;
		zrot << -1, 0, 0,
				0, -1, 0,
				0, 0, 1;
	}
};

static RotationMatrices rotators;

//always sort coordinates by mapnum
//require that all atoms of the core have a mapnum
//heavy atoms of the core only - full molecule is needed in case core has < 3 atoms
void ScaffoldIndexer::createCanonicalCoords(const Conformer& c, const Reaction::Decomposition& decomp,
		ECoords& coords, Orienter& orient) const
{
	ROMol& fullmol = c.getOwningMol();

	vector<MapNumInfo> connecting;
	vector<MapNumInfo> remaining;
	for (unsigned i = 0, n = decomp.core.size(); i < n; i++)
	{
		int idx = decomp.core[i];
		Atom *a = fullmol.getAtomWithIdx(idx);
		assert(a->hasProp(ATOM_MAP_NUM));
		assert(a->getAtomicNum() != 1);
		int mapnum = 0;
		a->getProp(ATOM_MAP_NUM, mapnum);

		if (connectingMapNums.count(mapnum) > 0)
			connecting.push_back(MapNumInfo(idx, mapnum));
		else
			remaining.push_back(MapNumInfo(idx, mapnum));
	}

	//if there are fewer than 3 total atoms, add first atom in connections
	if(connecting.size() + remaining.size() < 3)
	{
		assert(decomp.connections.size() > 0);
		assert(decomp.connections[0].size() > 0);
		const Reaction::Connection& c = decomp.connections[0][0];
		remaining.push_back(MapNumInfo(c.reactantIndex, c.reactantMap));
	}

	sort(connecting.begin(), connecting.end());
	sort(remaining.begin(), remaining.end());

	coords = ECoords::Zero(connecting.size() + remaining.size(), 3);

	//connecting always go first
	unsigned nc = connecting.size();
	for (unsigned i = 0; i < nc; i++)
	{
		RDGeom::Point3D pt = c.getAtomPos(connecting[i].idx);
		coords(i, 0) = pt.x;
		coords(i, 1) = pt.y;
		coords(i, 2) = pt.z;
	}
	for (unsigned i = 0, n = remaining.size(); i < n; i++)
	{
		RDGeom::Point3D pt = c.getAtomPos(remaining[i].idx);
		coords(i + nc, 0) = pt.x;
		coords(i + nc, 1) = pt.y;
		coords(i + nc, 2) = pt.z;
	}

	//center coordinates
	double nr = coords.rows();
	Vector3d center = coords.colwise().sum() / nr;
	coords.rowwise() -= center.transpose();
	orient.addTranslation(-center);

	//align to moment of inertia
	double Ixx = coords.block(0, 1, nr, 2).squaredNorm();
	double Iyy = coords.col(0).squaredNorm() + coords.col(2).squaredNorm();
	double Izz = coords.block(0, 0, nr, 2).squaredNorm();
	double Iyx = -coords.col(0).dot(coords.col(1));
	double Iyz = -coords.col(1).dot(coords.col(2));
	double Ixz = -coords.col(0).dot(coords.col(2));

	EMatrix3 I;
	I << Ixx, Iyx, Ixz,
			Iyx, Iyy, Iyz,
			Ixz, Iyz, Izz;

	SelfAdjointEigenSolver<EMatrix3> es(I);
	EMatrix3 principalAxes = es.eigenvectors();
	if (principalAxes.determinant() < 0) //has reflection
		principalAxes = -principalAxes;

	//rotation to align to prinicpal axes
	coords = coords * principalAxes;
	orient.addRotation(principalAxes);
	//next standardize alignment of moments

	if (!isPositiveBiased(coords, 1))
	{
		//rotate around x axis
		coords = coords * rotators.xrot;
		orient.addRotation(rotators.xrot);
	}
	if (!isPositiveBiased(coords, 0))
	{
		//rotate around y
		coords = coords * rotators.yrot;
		orient.addRotation(rotators.yrot);
	}
}

//for identifying best matches
struct RMSDSortInfo
{
	double rmsd;
	unsigned index;

	RMSDSortInfo() :
			rmsd(HUGE_VAL), index(0)
	{
	}
	RMSDSortInfo(double r, unsigned i) :
			rmsd(r), index(i)
	{
	}
	bool operator<(const RMSDSortInfo& rhs) const
	{
		return rmsd < rhs.rmsd;
	}
	bool operator==(const RMSDSortInfo& rhs) const
	{
		return rmsd == rhs.rmsd;
	}
};

//finds the best scaffold cluster for the passed coordinates of the core scaffold
//the cluster index is put in idx, with the closest match first
//if the best match does not meet the matching criteria, return false
bool ScaffoldIndexer::findBest(const ECoords& coords,
		vector<unsigned>& idx) const
		{
	idx.clear();

	//keep track of all valid matches
	vector<RMSDSortInfo> valid;
	//and the single very best match
	unsigned besti = 0;
	double bestrmsd = HUGE_VAL;
	for (unsigned i = 0, n = clusters.size(); i < n; i++)
	{
		double connectR = 0, totalR = 0;
		if (calcRMSDSquares(coords, clusters[i].center, connectR, totalR))
		{
			valid.push_back(RMSDSortInfo(totalR, i));
		}
		if (totalR < bestrmsd)
		{
			bestrmsd = totalR;
			besti = i;
		}
	}

	if (valid.size() > 0) //found something
	{
		sort(valid.begin(), valid.end());
		for (unsigned i = 0, n = valid.size(); i < n; i++)
		{
			idx.push_back(valid[i].index);
		}
		return true;
	}
	else
	{
		idx.push_back(besti);
		return false;
	}
}

void ScaffoldIndexer::addRotation(unsigned s, const ECoords& coords, Orienter& orient)
{
	EMatrix3 rot = computeRotation(clusters[s].center, coords);
	orient.addRotation(rot);
}


//add a new (unique) scaffold conformation
unsigned ScaffoldIndexer::addScaffold(const Conformer& c, const Reaction::Decomposition& decomp, Orienter& orient)
{
	vector<unsigned> idx;
	ECoords coords;
	createCanonicalCoords(core, coords, orient);
	if (!findBest(coords, idx))
	{
		//must create new scaffold
		clusters.push_back(ScaffoldInfo());
		ScaffoldInfo& info = clusters.back();
		info.center = coords;
		info.count++;
		return clusters.size() - 1;
	}
	else
	{
		assert(idx.size() > 0);
		unsigned index = idx[0];
		clusters[index].count++;
		addRotation(index, coords, orient);
		return index;
	}
}

void ScaffoldIndexer::dumpCounts(ostream& out) const
		{
	for (unsigned i = 0, n = clusters.size(); i < n; i++)
	{
		out << i << " " << clusters[i].count << "\n";
	}
}
