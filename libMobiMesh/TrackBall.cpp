//
//  TrackBall.cpp
//  MobiMesh
//
//  Created by Ong Yuh Shin on 4/29/11.
//  Copyright 2011 National University of Singapore. All rights reserved.
//

#include <iostream>

#include "glm.hpp"
#include "matrix_transform.hpp"

#include "TrackBall.h"

namespace MobiMesh {

bool TrackBall::get_rotation(float screen_x_init, float screen_y_init,
                             float screen_x_final, float screen_y_final,
                             float width, float height,
                             float& angle_radians,
                             float& axis_x, float& axis_y, float& axis_z)
{
    glm::vec3 initPos = map_to_sphere(screen_x_init, screen_y_init, width, height);
    glm::vec3 finalPos = map_to_sphere(screen_x_final, screen_y_final, width, height);

    glm::vec3 rotation_axis = glm::cross(initPos, finalPos);

    if (rotation_axis == glm::vec3(0.0))
        return false;

    axis_x = rotation_axis[0];
    axis_y = rotation_axis[1];
    axis_z = rotation_axis[2];

    angle_radians = std::acos(glm::dot(initPos, finalPos));

    return true;
}

glm::vec3 TrackBall::map_to_sphere(float screen_x, float screen_y,
                                   float width, float height)
{
    glm::vec3 sphere_coords;

	// project x, y onto a hemi-sphere centered at (0,0,0) with radius 1
    // first scale the viewport to [0,0]-[2,2].
	sphere_coords[0] = (2.0 * screen_x) / width;
	sphere_coords[1] = (2.0 * screen_y) / height;
    // translate (0,0) of viewport to center (0,0,0)
    sphere_coords[0] = sphere_coords[0] - 1;
    sphere_coords[1] = sphere_coords[1] - 1;

    // compute z based on radius of sphere (1.0)
    double z2 = 1 - sphere_coords[0]*sphere_coords[0] - sphere_coords[1]*sphere_coords[1];
    sphere_coords[2] = z2 > 0 ? sqrt(z2) : 0;

    // normalize the vector
    return glm::normalize(sphere_coords);
}

}