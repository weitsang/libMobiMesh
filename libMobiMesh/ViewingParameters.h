/*
 *  ViewingParameters.h
 *  MobiMesh
 *
 *  Created by Ong Yuh Shin on 3/10/11.
 *  Copyright 2011 National University of Singapore. All rights reserved.
 *
 */

#ifndef __VIEWING_PARAMETERS_H_
#define __VIEWING_PARAMETERS_H_

#include "VectorT.h"
#include "Plane3d.h"

namespace MobiMesh {

class ViewingParameters {

private:
	float           modelview_matrix_[16];
	float           fovy_,aspect_,tolerance_square_;
	OpenMesh::Vec3f eye_pos_,right_dir_,up_dir_,view_dir_;
	Plane3d         frustum_plane_[4];

public:

	ViewingParameters() { fovy_ = 45.0f; aspect_ = 1.0f; tolerance_square_ = 0.001f; }

	void increase_tolerance()                           { tolerance_square_ *= 5.0f; }
	void decrease_tolerance()                           { tolerance_square_ /= 5.0f; }  
	void set_fovy(float _fovy)                          { fovy_              = _fovy; }
	void set_aspect(float _aspect)                      { aspect_            = _aspect; }
	void set_tolerance_square(float _tolerance_square)  { tolerance_square_  = _tolerance_square; }

	// Note! We inline these get functions to speed up view-dependent queries; (0.5.7)

	inline float fovy()             const  { return  fovy_; }
	inline float aspect()           const  { return  aspect_; }
	inline float tolerance_square() const  { return  tolerance_square_; } 
	inline const OpenMesh::Vec3f& eye_pos()   const  { return  eye_pos_; }
	inline const OpenMesh::Vec3f& right_dir() const  { return  right_dir_; }
	inline const OpenMesh::Vec3f& up_dir()    const  { return  up_dir_; }
	inline const OpenMesh::Vec3f& view_dir()  const  { return  view_dir_; }
	inline OpenMesh::Vec3f& eye_pos()                { return  eye_pos_; }
	inline OpenMesh::Vec3f& right_dir()              { return  right_dir_; }
	inline OpenMesh::Vec3f& up_dir()                 { return  up_dir_; }
	inline OpenMesh::Vec3f& view_dir()               { return  view_dir_; }

	void frustum_planes( Plane3d _plane[4] ) {
		for (unsigned int i=0; i<4 ; ++i) _plane[i] = frustum_plane_[i];
	}
	void get_modelview_matrix(float _modelview_matrix[16]) {
		for (unsigned int i=0; i<16; ++i) _modelview_matrix [i] =  modelview_matrix_[i];
	}
	void set_modelview_matrix(const float _modelview_matrix[16]) {
		for (unsigned int i=0; i<16; ++i)  modelview_matrix_[i] = _modelview_matrix [i];   
	}

	// Note! We allow direct access to local ModelView Matrix; (0.5.6)
	float* modelview_matrix() { return &modelview_matrix_[0]; }

	// Note! We move this criteria into 'Class ViewingParameters' to avoid frequent memory copy; (0.5.6)
	/**
	 * Checks if a point is outside the view frustum.
	 *
	 * @param[in] pos coordinate of point to be checked for containment
	 * within view frustum
	 * @param[in] radius radius of bounding sphere for the point pos
	 * @return true iff the bounding sphere of pos falls at least
	 * partially within the view frustum
	 * @note view frustum is a 4-sided semi-infinite pyramid (not 6-sided)
	 */
	bool outside_view_frustum(const OpenMesh::Vec3f &pos, float radius) {  
		for (int i = 0; i < 4; ++i) { if (frustum_plane_[i].signed_distance(pos) < -radius) return true; } return false;
	}
	// Note! We move this function inside to be compatible with later migrations to iOS; (0.5.6)
	void update_viewing_configurations()
	{
		// |a11 a12 a13|-1       |  a33a22-a32a23  -(a33a12-a32a13)   a23a12-a22a13 |
		// |a21 a22 a23| = 1/DET*|-(a33a21-a31a23)   a33a11-a31a13  -(a23a11-a21a13)|
		// |a31 a32 a33|         |  a32a21-a31a22  -(a32a11-a31a12)   a22a11-a21a12 |
		//  DET  =  a11(a33a22-a32a23)-a21(a33a12-a32a13)+a31(a23a12-a22a13)

		float           invdet;  
		float           a11, a12, a13, a21, a22, a23, a31, a32, a33;
		OpenMesh::Vec3f inv_rot[3], trans;

		a11      = modelview_matrix_[0];  a21      = modelview_matrix_[1];  a31      = modelview_matrix_[2];
		a12      = modelview_matrix_[4];  a22      = modelview_matrix_[5];  a32      = modelview_matrix_[6];
		a13      = modelview_matrix_[8];  a23      = modelview_matrix_[9];  a33      = modelview_matrix_[10];
		trans[0] = modelview_matrix_[12]; trans[1] = modelview_matrix_[13]; trans[2] = modelview_matrix_[14];

		invdet=a11*(a33*a22-a32*a23) - a21*(a33*a12-a32*a13) + a31*(a23*a12-a22*a13);
		invdet= (float) 1.0/invdet;

		(inv_rot[0])[0] =  (a33*a22-a32*a23) * invdet; (inv_rot[0])[1] = -(a33*a12-a32*a13) * invdet; (inv_rot[0])[2] =  (a23*a12-a22*a13) * invdet;
		(inv_rot[1])[0] = -(a33*a21-a31*a23) * invdet; (inv_rot[1])[1] =  (a33*a11-a31*a13) * invdet; (inv_rot[1])[2] = -(a23*a11-a21*a13) * invdet;
		(inv_rot[2])[0] =  (a32*a21-a31*a22) * invdet; (inv_rot[2])[1] = -(a32*a11-a31*a12) * invdet; (inv_rot[2])[2] =  (a22*a11-a21*a12) * invdet;

		eye_pos_   = - OpenMesh::Vec3f(dot(inv_rot[0], trans), dot(inv_rot[1], trans), dot(inv_rot[2], trans));
		right_dir_ =   OpenMesh::Vec3f(a11, a12, a13);
		up_dir_    =   OpenMesh::Vec3f(a21, a22, a23);
		view_dir_  = - OpenMesh::Vec3f(a31, a32, a33);

		OpenMesh::Vec3f normal[4];  
		//float aspect = width() / height();
		float half_theta = fovy() * 0.5f;
		float half_phi = atanf(aspect() * tanf(half_theta));

		float sin1 = sinf(half_theta), cos1 = cosf(half_theta);
		float sin2 = sinf(half_phi),   cos2 = cosf(half_phi);

		normal[0] =  cos2 * right_dir_ + sin2 * view_dir_;
		normal[1] = -cos1 * up_dir_    - sin1 * view_dir_;
		normal[2] = -cos2 * right_dir_ + sin2 * view_dir_;
		normal[3] =  cos1 * up_dir_    - sin1 * view_dir_;

		for (int i=0; i<4; i++) frustum_plane_[i] = Plane3d(normal[i], eye_pos_);
	}		
};

}

#endif