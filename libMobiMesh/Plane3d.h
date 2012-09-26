/*
 *  Plane3d.h
 *  MobiMesh
 *
 *  Created by Ong Yuh Shin on 3/10/11.
 *  Copyright 2011 National University of Singapore. All rights reserved.
 *
 */

#ifndef __PLANE3D_H_
#define __PLANE3D_H_

#include "VectorT.h"

namespace MobiMesh {

/**
 * Represents a plane in 3D space with some equation of form
 * ax + by + cz + d = 0
 * where (a,b,c) is the normal vector, (x,y,z) is a point on the plane.
 *
 * @TODO make n_ and d_ private and provide suitable accessors.
 */
class Plane3d
{

public:
    typedef OpenMesh::Vec3f          vector_type;
    typedef vector_type::value_type  value_type;

	Plane3d() : d_(0) { }
	Plane3d(const vector_type &_dir, const vector_type &_pnt)
		: n_(_dir), d_(0) { n_.normalize(); d_ = -dot(n_,_pnt); }

	value_type signed_distance(const OpenMesh::Vec3f &_p)
	{
		return  dot(n_ , _p) + d_;
	}

    vector_type n_;
    value_type  d_;
};

}

#endif