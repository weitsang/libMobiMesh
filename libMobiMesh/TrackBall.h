//
//  TrackBall.h
//  MobiMesh
//
//  Created by Ong Yuh Shin on 4/29/11.
//  Copyright 2011 National University of Singapore. All rights reserved.
//

#ifndef __TRACK_BALL_H_
#define __TRACK_BALL_H_

#include "glm.hpp"

namespace MobiMesh {

/**
 * Class to perform TrackBall rotation.
 * @see http://viewport3d.com/trackball.htm
 */
class TrackBall
{
public:

    /**
     * Create the trackball object to compute the rotation angle and
     * axis in 3D object space when the touch/mouse moves from an initial
     * to a final point on the viewport.
     */
    TrackBall() { }

    /**
     * Get the rotation angle and axis when the mouse/touch is
     * moved from one point to another. Note that the bottom left corner
     * of the screen is (0,0) and the x value increases as we move right,
     * the y value increases as we move up
     *
     * @param[in] screen_x_init initial x coordinate on screen/viewport
     * @param[in] screen_y_init initial y coordinate on screen/viewport
     * @param[in] screen_x_final final x coordinate on screen/viewport
     * @param[in] screen_y_final final y coordinate on screen/viewport
     * @param[in] width viewport width
     * @param[in] height viewport height
     * @param[out] angle_radians computed angle of rotation
     * @param[out] axis_x x coordinate of the computed axis of rotation
     * @param[out] axis_y y coordinate of the computed axis of rotation
     * @param[out] axis_z z coordinate of the computed axis of rotation
     * @return true iff the input results in a valid rotation
     */
    bool get_rotation(float screen_x_init, float screen_y_init,
                      float screen_x_final, float screen_y_final,
                      float width, float height,
                      float& angle_radians,
                      float& axis_x, float& axis_y, float& axis_z);

private:

    TrackBall(const TrackBall&);
    TrackBall& operator=(const TrackBall&);

    /**
     * Maps the screen coordinates to a unit sphere centered at (0,0,0)
     */
    glm::vec3 map_to_sphere(float screen_x, float screen_y,
                            float width, float height);
};

}

#endif