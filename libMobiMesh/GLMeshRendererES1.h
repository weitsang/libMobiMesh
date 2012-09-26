/*
 *  GLMeshRenderer.h
 *  MobiMesh
 *
 *  Created by Ong Yuh Shin on 3/11/11.
 *  Copyright 2011 National University of Singapore. All rights reserved.
 *
 */

#ifndef __GL_MESH_RENDERER_ES1_H_
#define __GL_MESH_RENDERER_ES1_H_

#include <vector>

#include <OpenGLES/ES1/gl.h>
#include <OpenGLES/ES1/glext.h>
#include <OpenGLES/ES2/gl.h>

#include "IMesh.h"
#include "IMeshRenderer.h"
#include "GLMeshRenderer.h"
#include "GLMeshRendererViewingParameters.h"

namespace MobiMesh {

/**
 * Concrete class responsible for making OpenGL ES 1.1 specific calls
 * when rendering.
 */
class GLMeshRendererES1 : public GLMeshRenderer
{
public:

	///////////// Constructors

    /**
     * Creates a renderer that uses OpenGL ES 1.1 to render
     * @param[in] width the viewport width
     * @param[in] height the viewport height
     */
	GLMeshRendererES1(GLsizei width, GLsizei height);

    /**
     * Creates a renderer that uses OpenGL ES 1.1 to render the supplied mesh
     * @param[in] mesh mesh object to be rendered. Caller is responsible for
     * ensuring that the lifetime of this object exceeds that of the created
     * renderer.
     * @param[in] width the viewport width
     * @param[in] height the viewport height
     */
	GLMeshRendererES1(IMesh* mesh, GLsizei width, GLsizei height);

	virtual ~GLMeshRendererES1();

protected:

    ///////////// Methods inherited from GLMeshRenderer

    void setup_modelview_projection();
    void setup_vertices(const GLvoid* pointer);
    void setup_normals(const GLvoid* pointer);
    void setup_colors(const GLvoid* pointer);
    void rendering_setup();
    void rendering_cleanup();

private:

    GLMeshRendererES1(const GLMeshRendererES1&);
    GLMeshRendererES1& operator=(const GLMeshRendererES1&);
};

}

#endif
