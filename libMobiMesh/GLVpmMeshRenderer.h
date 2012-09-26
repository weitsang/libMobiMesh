/*
 *  GLVpmMeshRenderer.h
 *  MobiMesh
 *
 *  Created by Ong Yuh Shin on 3/29/11.
 *  Copyright 2011 National University of Singapore. All rights reserved.
 *
 */

#ifndef __GL_VPM_MESH_RENDERER_H_
#define __GL_VPM_MESH_RENDERER_H_

#include "GLViewDependentMeshRenderer.h"
#include "ViewDependentProgressiveMesh.h"

namespace MobiMesh {

/**
 * Concrete class for adapting and rendering .vpm meshes.
 */
class GLVpmMeshRenderer : public GLViewDependentMeshRenderer
{
public:

	///////////// Constructors

	/**
     * Creates a renderer that refines the supplied VPM mesh based on the
     * modelview matrix.
     * @param[in] renderer renderer object to use for rendering the mesh using
     * the appropriate OpenGL ES version. Caller must ensure that lifetime of
     * this object outlives that of the created GLVpmMeshRenderer object.
     * @param[in] mesh the mesh to be refined and rendered. This mesh will
     * replace whatever mesh is currently used by the supplied renderer (if
     * any). Thus, if necessary, the caller should free up memory used by
     * renderer's current mesh.
	 */
	GLVpmMeshRenderer(GLMeshRenderer& renderer, ViewDependentProgressiveMesh* mesh);

	~GLVpmMeshRenderer();

    ///////////// Methods inherited from GLViewDependentMeshRenderer

    void activate();
    void deactivate();

protected:

    bool adapt_mesh();

private:

	GLVpmMeshRenderer(const GLVpmMeshRenderer&);
	GLVpmMeshRenderer& operator=(const GLVpmMeshRenderer&);

	ViewDependentProgressiveMesh* mesh_;
    ViewingParameters viewing_params_;

    bool need_adapt_;
};

}

#endif
