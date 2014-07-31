/*
 * Orienter.h
 *
 *  Created on: Jul 29, 2014
 *      Author: dkoes
 *
 *  Maintains a translation vector and rotation matrix for reorienting coordinates.
 */

#ifndef ORIENTER_H_
#define ORIENTER_H_

#include <eigen3/Eigen/Core>
#include <vector>
#include <rdkit/Geometry/point.h>

typedef Eigen::Matrix<double, Eigen::Dynamic, 3> ECoords;
typedef Eigen::Matrix3d EMatrix3; //3x3 double
typedef Eigen::Vector3d EVector3;

class Orienter
{
	EVector3 translate; //gets added to coords
	EMatrix3 rotate; //then multiply by this - actually transpose of rotation matrix for right multiply
	public:
	Orienter() :
			translate(Eigen::Vector3d::Zero()), rotate(
					Eigen::Matrix3d::Identity())
	{

	}

	//accumlate a translation vector
	void addTranslation(const EVector3& t)
	{
		translate += t;
	}

	void addRotation(const EMatrix3& r)
	{
		rotate *= r;
	}

	void reorient(ECoords& coords) const
	{
		coords.rowwise() += translate.transpose();
		coords *= rotate;
	}

	//reorient rd points
	void reorient(std::vector<RDGeom::Point3D>& coords) const
	{
		for(unsigned i = 0, n = coords.size(); i < n; i++)
		{
			EVector3 pt = EVector3(coords[i].x, coords[i].y, coords[i].z);
			pt = (pt + translate).transpose()*rotate;
			coords[i].x = pt(0);
			coords[i].y = pt(1);
			coords[i].z = pt(2);
		}
	}
};

#endif /* ORIENTER_H_ */
