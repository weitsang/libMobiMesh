//
//  GLMeshRendererViewingParameters.cpp
//  MobiMesh
//
//  Created by Ong Yuh Shin on 4/10/11.
//  Copyright 2011 National University of Singapore. All rights reserved.
//

#include "GLMeshRendererViewingParameters.h"

#include "glm.hpp"
#include "matrix_transform.hpp"

namespace MobiMesh {

GLMeshRendererViewingParameters::GLMeshRendererViewingParameters()
	: eye_x_(0.0), eye_y_(0.0), eye_z_(0.0),
    lookat_x_(0.0), lookat_y_(0.0), lookat_z_(0.0),
    up_x_(0.0), up_y_(0.0), up_z_(0.0),
    bounding_length_(0.0), min_distance_(0.0), max_distance_(0.0),
    center_x_(0.0), center_y_(0.0), center_z_(0.0),
    mesh_(NULL)
{
    set_default_view();
}

GLMeshRendererViewingParameters::GLMeshRendererViewingParameters(
    IMesh* mesh)
    : eye_x_(0.0), eye_y_(0.0), eye_z_(0.0),
    lookat_x_(0.0), lookat_y_(0.0), lookat_z_(0.0),
    up_x_(0.0), up_y_(0.0), up_z_(0.0),
    bounding_length_(0.0), min_distance_(0.0), max_distance_(0.0),
    center_x_(0.0), center_y_(0.0), center_z_(0.0),
    mesh_(mesh)
{
    set_default_view();
}

void GLMeshRendererViewingParameters::set_default_view()
{
    recompute_bounding_box();

	eye_x_ = center_x_;
	eye_y_ = center_y_;
	eye_z_ = center_z_ + get_eye_to_lookat_distance();

    lookat_x_ = center_x_;
    lookat_y_ = center_y_;
    lookat_z_ = center_z_;

    up_x_ = 0.0;
    up_y_ = 1.0;
    up_z_ = 0.0;

	min_distance_    = 0.01 * get_eye_to_lookat_distance();
	max_distance_    =  100 * get_eye_to_lookat_distance();

    lookat_matrix_ = glm::lookAt(
        glm::vec3(eye_x_, eye_y_, eye_z_),
        glm::vec3(lookat_x_, lookat_y_, lookat_z_),
        glm::vec3(up_x_, up_y_, up_z_));

    transformation_matrix_ = glm::mat4(1.0);
}

void GLMeshRendererViewingParameters::recompute_bounding_box()
{
    // Compute bounding box (of mesh/vertices)
    float max_length;
    float min_x = 0,    min_y = 0,    min_z = 0;
    float max_x = 0,    max_y = 0,    max_z = 0;

    if (!mesh_ ||
        (!mesh_->has_face_vertices() && !mesh_->has_vertex_indices()))
    {
        center_x_ = center_y_ = center_z_ = 0;
        bounding_length_ = 0;
        return;
    }

    const IMesh::PointContainer* points = mesh_->points();
    if (!points)
        return;

    if (points->size() == 0)
        return;

    center_x_ = min_x = max_x = (*points)[0][0];
    center_y_ = min_y = max_y = (*points)[0][1];
    center_z_ = min_z = max_z = (*points)[0][2];

    for (size_t vi = 1; vi < points->size(); ++vi)
    {
        float x = (*points)[vi][0];
        float y = (*points)[vi][1];
        float z = (*points)[vi][2];
        min_x = min_x > x ? x : min_x;
        min_y = min_y > y ? y : min_y;
        min_z = min_z > z ? z : min_z;
        max_x = max_x < x ? x : max_x;
        max_y = max_y < y ? y : max_y;
        max_z = max_z < z ? z : max_z;
    }
    center_x_ = (max_x + min_x)/2;
    center_y_ = (max_y + min_y)/2;
    center_z_ = (max_z + min_z)/2;
    max_length = max_x - min_x;

    if (max_y - min_y > max_length)
        max_length = max_y - min_y;
    if (max_z - min_z > max_length)
        max_length = max_z - min_z;

    bounding_length_ = max_length;
}

void GLMeshRendererViewingParameters::translate(float dx, float dy, float dz)
{
    glm::mat4 temp_matrix = glm::translate(glm::mat4(1.0),
                                           glm::vec3(dx, dy, dz));
    transformation_matrix_ = temp_matrix * transformation_matrix_;
}

void GLMeshRendererViewingParameters::rotate(
    float angle_degrees, float axis_x, float axis_y, float axis_z)
{
    glm::vec3 t = glm::vec3(center_x_, center_y_, center_z_);
    glm::mat4 temp_matrix = glm::translate(glm::mat4(1.0), t);
    temp_matrix = glm::rotate(temp_matrix,
                              angle_degrees,
                              glm::vec3(axis_x, axis_y, axis_z));
    temp_matrix = glm::translate(temp_matrix, -t);
    transformation_matrix_ = temp_matrix * transformation_matrix_;
}

void GLMeshRendererViewingParameters::scale(float scale_factor)
{
    glm::vec3 t = glm::vec3(center_x_, center_y_, center_z_);
    glm::mat4 temp_matrix = glm::translate(glm::mat4(1.0), t);
    temp_matrix = glm::scale(temp_matrix,
                             glm::vec3(scale_factor, scale_factor, scale_factor));
    temp_matrix = glm::translate(temp_matrix, -t);
    transformation_matrix_ = temp_matrix * transformation_matrix_;
}

}