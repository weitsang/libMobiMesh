/*
 *  IMeshRenderer.h
 *  MobiMesh
 *
 *  Created by Ong Yuh Shin on 3/8/11.
 *  Copyright 2011 National University of Singapore. All rights reserved.
 *
 */

#ifndef __IMESH_RENDERER_H_
#define __IMESH_RENDERER_H_

#include <vector>

#include "VectorT.h"
#include "glm.hpp"

namespace MobiMesh {

class IMesh;

/**
 * Renderer interface: to render a single mesh object.
 * It supports the usual rotate, scale and translate methods.
 * It is also responsible for setting up a view that places the
 * mesh in the center of the screen.
 * The modelview matrix is made up of a lookat * transformation
 * matrix, where the lookat matrix represents the position of the
 * camera, and the transformation matrix is a cumulative transformation
 * on the mesh being rendered. (All transformations are cumulative,
 * i.e. performed on the existing transformation matrix)
 */
class IMeshRenderer
{
public:

	///////////// Constructors

	virtual ~IMeshRenderer() { }

	///////////// Rendering
	
	/**
	 * Renders the mesh. The implementing classes should store any
	 * necessary rendering parameters as its data members.
	 */
	virtual void render() =0;

	virtual void translate(float dx, float dy, float dz) =0;

	virtual void rotate(float angle_degrees,
                        float axis_x, float axis_y, float axis_z) =0;

    virtual void scale(float scale_factor) =0;

    virtual void set_default_view() =0;

    // Modelview matrix = lookat * transformation
    virtual void capture_view(glm::mat4& lookat,
                              glm::mat4& transformation) const =0;
    virtual void set_view(const glm::mat4& lookat,
                          const glm::mat4 transformation) =0;

    /// Longest length of bounding box in object space
    /// (estimate of the size of mesh in object space)
    virtual float get_bounding_length() const =0;

    virtual int get_viewport_width() const =0;
    virtual int get_viewport_height() const =0;

    /**
     * Activate the renderer. This sets up the necessary resources
     * for rendering the mesh on screen (e.g. buffers, mesh, etc.)
     */
    virtual void activate() =0;

    /**
     * Deactivates the renderer. This frees up resources that are
     * not required for rendering on screen. It retains what is
     * necessary for activate() to reload the mesh quickly.
     */
    virtual void deactivate() =0;

    /**
     * Renderers should default to active upon construction.
     * All calls in this interface besides activate() have undefined
     * behavior if the renderer is not active.
     * @return true iff renderer is active
     */
    virtual bool is_active() =0;
};

}

#endif
