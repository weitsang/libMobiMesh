//
//  GLTvpmMeshRenderer.h
//  MobiMesh
//
//  Created by Ong Yuh Shin on 5/19/11.
//  Copyright 2011 National University of Singapore. All rights reserved.
//

#ifndef __GL_TVPM_MESH_RENDERER_H_
#define __GL_TVPM_MESH_RENDERER_H_

#include "GLViewDependentMeshRenderer.h"
#include "TrulyViewDependentPM.h"

namespace MobiMesh {

/**
 * Concrete class for rendering and adapting .tvpm and .tvpm64 meshes.
 */
class GLTvpmMeshRenderer : public GLViewDependentMeshRenderer
{
public:

	///////////// Constructors

	/**
     * Creates a renderer that refines the supplied TVPM mesh based on the
     * percentage of the viewport the mesh's visible faces occupy.
     *
     * @param[in] renderer renderer object to use for rendering the mesh
     * using the appropriate OpenGL ES version. Caller must ensure that
     * lifetime of
     * this object outlives that of the created GLTvpmMeshRenderer object.
     * @param[in] mesh the mesh to be refined and rendered. This mesh will
     * replace whatever mesh is currently used by the supplied renderer (if
     * any). Thus, if necessary, the caller should free up memory used by
     * renderer's current mesh.
	 */
	GLTvpmMeshRenderer(GLMeshRenderer& renderer,
                       TrulyViewDependentPM* mesh);

	~GLTvpmMeshRenderer();

    ///////////// Methods inherited from GLViewDependentMeshRenderer
    
    void activate();
    void deactivate();

protected:

    bool adapt_mesh();

private:
    
	GLTvpmMeshRenderer(const GLTvpmMeshRenderer&);
	GLTvpmMeshRenderer& operator=(const GLTvpmMeshRenderer&);

    void initialize();

    /**
     * Render the mesh with each color assigned to each face and reads
     * the pixels back to determine which faces are visible.
     * The value read is stored in the image_ array
     */
    void read_pixels();

    /**
     * Based on the values in image_ array, this method computes and
     * updates the vertex weights in the tvpm mesh_.
     */
    void update_vertex_weights();

	TrulyViewDependentPM* mesh_;

    // To store pixel value of mesh
    unsigned char* image_;

    bool read_pixels_, need_further_refine_;

    // Number of times the pixels have been read for the current view
    int read_pixels_count_;
    // Edge collapse threshold (constant wrt to viewport size)
    int ec_threshold_;

    // Number of state of the renderer.
    static const int NUM_STATES = 3;

    // Describes the vertex threshold to be used for each state.
    // This must be initialized with the vertex split threshold for each
    // possible state. The vertex split threshold is given in terms of
    // a percentage of the viewport area (rather than a specific weight)
    // Values are ordered from large to small.
    static const float vertex_threshold[NUM_STATES];

    // Current state that the renderer is in, needs to be set before
    // adapting the mesh so that the appropriate vertex threshold can
    // be passed. State should be an integer referring to something in State
    int curr_state;
};

}

#endif