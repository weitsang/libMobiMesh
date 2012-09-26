/*
 *  GLMeshRenderer.cpp
 *  MobiMesh
 *
 *  Created by Ong Yuh Shin on 3/11/11.
 *  Copyright 2011 National University of Singapore. All rights reserved.
 *
 */

#include <algorithm>

#include "glm.hpp"
#include "matrix_transform.hpp"

#include "GLMeshRendererES1.h"

namespace MobiMesh {

GLMeshRendererES1::GLMeshRendererES1(
    GLsizei viewport_width, GLsizei viewport_height)
    : GLMeshRenderer(viewport_width, viewport_height)
{ }

GLMeshRendererES1::GLMeshRendererES1(
    IMesh*  mesh,
    GLsizei viewport_width,
    GLsizei viewport_height)
	: GLMeshRenderer(mesh, viewport_width, viewport_height)
{ }

GLMeshRendererES1::~GLMeshRendererES1()
{ }

    
void GLMeshRendererES1::setup_modelview_projection()
{
    // Setup the projection matrix
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glViewport(0, 0, get_viewport_width(), get_viewport_height());

    glm::mat4 perspective_matrix =
        glm::perspective((GLfloat) GLMeshRenderer::FOVY,
                         (GLfloat) get_viewport_width() / (GLfloat) get_viewport_height(),
                         view_params_.get_near_plane(),
                         view_params_.get_far_plane());

    glMultMatrixf(&perspective_matrix[0][0]);

    // Setup the modelview matrix
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    glm::mat4 lookat, transformation;
    capture_view(lookat, transformation); // get the current view
    glm::mat4 modelview_matrix = lookat * transformation;

    glMultMatrixf(&modelview_matrix[0][0]);
}

void GLMeshRendererES1::setup_vertices(const GLvoid *pointer)
{
    glVertexPointer(3, GL_FLOAT, 0, pointer);
}

void GLMeshRendererES1::setup_normals(const GLvoid *pointer)
{
    glNormalPointer(GL_FLOAT, 0, pointer);
}

void GLMeshRendererES1::setup_colors(const GLvoid *pointer)
{
    glColorPointer(4, GL_UNSIGNED_BYTE, 0, pointer);
}

void GLMeshRendererES1::rendering_setup()
{
    glEnableClientState(GL_VERTEX_ARRAY);
    if (render_with_color())
    {
        glDisable (GL_LIGHTING);
        glDisable (GL_DITHER);
        glEnableClientState(GL_COLOR_ARRAY);
    }
    else
    {
        glEnableClientState(GL_NORMAL_ARRAY);
    }
}

void GLMeshRendererES1::rendering_cleanup()
{
    glDisableClientState(GL_VERTEX_ARRAY);

    if (render_with_color())
    {
        glDisableClientState (GL_COLOR_ARRAY);
        glEnable(GL_DITHER);
        glEnable(GL_LIGHTING);
    }
    else
    {
        glDisableClientState(GL_NORMAL_ARRAY);
    }
}

}
