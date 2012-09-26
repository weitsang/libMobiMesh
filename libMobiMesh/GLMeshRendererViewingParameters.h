//
//  GLMeshRendererViewingParameters.h
//  MobiMesh
//
//  Created by Ong Yuh Shin on 4/10/11.
//  Copyright 2011 National University of Singapore. All rights reserved.
//

#ifndef __GL_MESH_RENDERER_VIEWING_PARAMETERS_H_
#define __GL_MESH_RENDERER_VIEWING_PARAMETERS_H_

#include "IMesh.h"
#include "glm.hpp"

namespace MobiMesh {

/**
 * Represents the viewing parameters when rendering with OpenGL.
 * Viewing parameters include the x,y panning, x,y rotation,
 * scale factor and also the mesh to be rendered.
 *
 * All transformations transformation calls are performed on
 * based on the stored modelview matrix.
 * The modelview matrix is madeup of lookat*transformation.
 * All transformations are pre-multiplied by the transformation matrix,
 * i.e. composed incrementally.
 * The lookat matrix only changes when the bounding box is changed.
 */
class GLMeshRendererViewingParameters
{
public:

    GLMeshRendererViewingParameters();

    /**
     * Initializes the set of viewing parameters with the given mesh
     * object. Caller is responsible to ensure that the lifetime of
     * the supplied mesh exceeds the lifetime of the
     * GLMeshRendererViewingParameters object being created.
     */
    explicit GLMeshRendererViewingParameters(IMesh* mesh);

    const IMesh* get_mesh() const { return mesh_; }
    void set_mesh(IMesh* mesh) { mesh_ = mesh; }

    void translate(float dx, float dy, float dz);

    void rotate(float angle_degrees,
                float axis_x, float axis_y, float axis_z);

    void scale(float scale_factor);

    float get_eye_x() const { return eye_x_; }
    float get_eye_y() const { return eye_y_; }
    float get_eye_z() const { return eye_z_; }

    float get_lookat_x() const { return lookat_x_; }
    float get_lookat_y() const { return lookat_y_; }
    float get_lookat_z() const { return lookat_z_; }

    // Up vector for the eye and lookat points (points along y-axis)
    float get_up_x() const { return up_x_; }
    float get_up_y() const { return up_y_; }
    float get_up_z() const { return up_z_; }

    float get_near_plane() const { return min_distance_; }
    float get_far_plane() const { return max_distance_; }

    /// Longest length of bounding box
    float get_bounding_length() const { return bounding_length_; }

	/**
	 * Recomputes viewing parameters to the default, and sets the modelview
     * matrix accordingly.
     * The bounding box is recomputed based on the current set of points.
     * Let (x, y, z) be the center of the bounding box, the default eye
     * position is (x, y, z + 1.5 * bounding_length_).
     * Let (eyex, eyey, eyez) be the eye position from before, the default
     * lookat position is (eyex, eyey, eyez - 2 * bounding_length_)
	 */
    void set_default_view();

    glm::mat4 get_modelview_matrix() const
    { return lookat_matrix_ * transformation_matrix_; }

    glm::mat4 get_lookat_matrix() const { return lookat_matrix_; }
    void set_lookat_matrix(const glm::mat4& lookat) { lookat_matrix_ = lookat; }

    glm::mat4 get_transformation_matrix() const { return transformation_matrix_; }
    void set_transformation_matrix(const glm::mat4& transformation)
    { transformation_matrix_ = transformation; }

private:

    // Mesh data used for rendering
    IMesh* mesh_;

    glm::mat4 transformation_matrix_;
    glm::mat4 lookat_matrix_;

    // Eye position
    float eye_x_, eye_y_, eye_z_;
    // Look at point
    float lookat_x_, lookat_y_, lookat_z_;
    // The current up direction
    float up_x_, up_y_, up_z_;

    // Longest length of bounding box along one of the 3 axes
    float bounding_length_;
    // Center of bounding box
    float center_x_, center_y_, center_z_;

    // Distance from eye of front and back of frustum plane
    // i.e. the near and far clipping planes
    float min_distance_, max_distance_;

    float get_eye_to_lookat_distance() const { return 1.5 * bounding_length_; }

    /**
     * Recompute the bounding box parameters based on the current set of points
     */
    void recompute_bounding_box();

    GLMeshRendererViewingParameters(const GLMeshRendererViewingParameters&);
    GLMeshRendererViewingParameters& operator=(const GLMeshRendererViewingParameters&);
};

}

#endif