//
//  GLMeshRenderer.cpp
//  MobiMesh
//
//  Created by Ong Yuh Shin on 5/5/11.
//  Copyright 2011 National University of Singapore. All rights reserved.
//

#include <algorithm>

#include "glm.hpp"
#include "matrix_transform.hpp"

#include "MobiMeshConstants.h"
#include "GLMeshRenderer.h"
#include "InvalidGeometryException.h"
#include "InactiveRendererException.h"

namespace MobiMesh {

GLMeshRenderer::GLMeshRenderer(GLsizei viewport_width,
                               GLsizei viewport_height)
    : with_color_(false),
    viewport_width_(viewport_width), viewport_height_(viewport_height),
    is_active_(false)
{
    activate();
}

GLMeshRenderer::GLMeshRenderer(IMesh*  mesh,
                               GLsizei viewport_width,
                               GLsizei viewport_height)
	: view_params_(mesh), with_color_(false),
    viewport_width_(viewport_width), viewport_height_(viewport_height),
    is_active_(false)
{
    activate();
}

GLMeshRenderer::~GLMeshRenderer()
{
    deactivate();
}

void GLMeshRenderer::set_color(const std::vector<OpenMesh::Vec4uc>& color)
{
#ifdef DEBUG
    if (!is_active())
    {
        throw InactiveRendererException("Renderer is inactive, call activate() "
                                        "to render");
    }
#endif

    glBindBuffer(GL_ARRAY_BUFFER, vbo_colors_);
    glBufferData(GL_ARRAY_BUFFER,
                 color.size() * sizeof(OpenMesh::Vec4uc),
                 &color[0],
                 GL_STATIC_DRAW);
}

void GLMeshRenderer::render()
{
#ifdef DEBUG
    if (!is_active())
    {
        throw InactiveRendererException("Renderer is inactive, call activate() "
                                        "to render");
    }
#endif

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glViewport(0, 0, get_viewport_width(), get_viewport_height());
    setup_modelview_projection();

    const IMesh* mesh = view_params_.get_mesh();

    if (!mesh) // no mesh to render
        return;
    else if (mesh->has_vertex_indices())
        bind_and_draw_vbo_draw_elements();
    else if (mesh->has_face_vertices())
        bind_and_draw_vbo_draw_arrays();
}

void GLMeshRenderer::activate()
{
    if (is_active())
        return;

    glGenBuffers(1, &vbo_vertices_);
    glGenBuffers(1, &vbo_normals_);
    glGenBuffers(1, &vbo_faces_);
    glGenBuffers(1, &vbo_colors_);

    is_active_ = true;

    // We can only update the mesh if the renderer is active so this
    // call comes after is_active_ is set to true. This will restore
    // the vertex/normal/face buffers to the state from before.
    // The color buffer is not restored. The caller has to restore it
    // by calling set_colors again if necessary.
    update_mesh();
}

void GLMeshRenderer::deactivate()
{
    if (!is_active())
        return;

    is_active_ = false;
    render_with_color(false);

    glDeleteBuffers(1, &vbo_vertices_);
    glDeleteBuffers(1, &vbo_normals_);
    glDeleteBuffers(1, &vbo_faces_);
    glDeleteBuffers(1, &vbo_colors_);
}

void GLMeshRenderer::update_mesh()
{
#ifdef DEBUG
    if (!is_active())
    {
        throw InactiveRendererException("Renderer is inactive, call activate() "
                                        "to render");
    }
#endif

    const IMesh* mesh = view_params_.get_mesh();

    if (!mesh) // no mesh to render
        return;
    else if (mesh->has_vertex_indices())
    {
        const IMesh::PointContainer* points =
            view_params_.get_mesh()->points();
        const IMesh::VertexNormalContainer* vertex_normals =
            view_params_.get_mesh()->vertex_normals();
        const IMesh::FaceContainer* faces =
            view_params_.get_mesh()->faces();

        if (!points || !vertex_normals || !faces)
            throw InvalidGeometryException("points/normals/faces do not exist");
        if (points->size() != vertex_normals->size())
        {
            throw InvalidGeometryException("number of vertex normals different "
                                           "from number of vertices and do not "
                                           "correspond");
        }

        num_faces_ = mesh->num_faces();

        glBindBuffer(GL_ARRAY_BUFFER, vbo_vertices_);
        glBufferData(GL_ARRAY_BUFFER,
                     points->size() * sizeof(OpenMesh::Vec3f),
                     &(*points)[0],
                     GL_STREAM_DRAW);

        glBindBuffer(GL_ARRAY_BUFFER, vbo_normals_);
        glBufferData(GL_ARRAY_BUFFER,
                     vertex_normals->size() * sizeof(OpenMesh::Vec3f),
                     &(*vertex_normals)[0],
                     GL_STREAM_DRAW);

        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, vbo_faces_);
        glBufferData(GL_ELEMENT_ARRAY_BUFFER,
                     faces->size() * sizeof(unsigned short),
                     &(*faces)[0],
                     GL_STREAM_DRAW);
    }
    else if (mesh->has_face_vertices())
    {
        const IMesh::PointContainer* vbo_points =
            view_params_.get_mesh()->face_vertices();
        const IMesh::VertexNormalContainer* vbo_normals =
            view_params_.get_mesh()->face_normals();

        if (!vbo_points || !vbo_normals)
            throw InvalidGeometryException("points/normals/faces do not exist");
        if (vbo_points->size() != vbo_normals->size())
        {
            throw InvalidGeometryException("number of vertex normals different "
                                           "from number of vertices and do not "
                                           "correspond");
        }

        num_faces_ = mesh->num_faces();

        glBindBuffer(GL_ARRAY_BUFFER, vbo_vertices_);
        glBufferData(GL_ARRAY_BUFFER,
                     3*num_faces_*sizeof(OpenMesh::Vec3f),
                     &(*vbo_points)[0],
                     GL_STREAM_DRAW);

        glBindBuffer(GL_ARRAY_BUFFER, vbo_normals_);
        glBufferData(GL_ARRAY_BUFFER,
                     3*num_faces_*sizeof(OpenMesh::Vec3f),
                     &(*vbo_normals)[0],
                     GL_STREAM_DRAW);
    }
}

void GLMeshRenderer::bind_and_draw_vbo_draw_elements()
{
    rendering_setup();

    glBindBuffer(GL_ARRAY_BUFFER, vbo_vertices_);
    setup_vertices(0);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, vbo_faces_);

    if (render_with_color())
    {
        glBindBuffer(GL_ARRAY_BUFFER, vbo_colors_);
        setup_colors(0);
    }
    else
    {
        glBindBuffer(GL_ARRAY_BUFFER, vbo_normals_);
        setup_normals(0);
    }

    glDrawElements(GL_TRIANGLES, num_faces_ * 3, GL_UNSIGNED_SHORT, 0);

    rendering_cleanup();
}

void GLMeshRenderer::render_with_color(bool color) 
{ 
    with_color_ = color;
    if (with_color_) {
        glBindBuffer(GL_ARRAY_BUFFER, vbo_colors_);
        setup_colors(0);
    } else {
        glBindBuffer(GL_ARRAY_BUFFER, vbo_normals_);
        setup_normals(0);
    }
}
    
void GLMeshRenderer::bind_and_draw_vbo_draw_arrays()
{
    rendering_setup();

    glBindBuffer(GL_ARRAY_BUFFER, vbo_vertices_);
    setup_vertices(0);

    // should not do this everytime.  Moved to bind/setup when
    // color changes - ooiwt
    //if (render_with_color())
    //{
    //    glBindBuffer(GL_ARRAY_BUFFER, vbo_colors_);
    //    setup_colors(0);
   // }
    //else
    //{
     //   glBindBuffer(GL_ARRAY_BUFFER, vbo_normals_);
    //    setup_normals(0);
    //}

    // Since each glDrawArrays call can only draw up to MAX_FACES*3 elements
    // each time, we split these draw calls into batches
    int start_index = 0;
    int num_draws = num_faces_ / Constants::MAX_FACES;
    int remainder = num_faces_ % Constants::MAX_FACES;
    for (int draw_count = 0; draw_count < num_draws; ++draw_count)
    {
        glDrawArrays(GL_TRIANGLES, start_index, Constants::MAX_FACES * 3);
        start_index += Constants::MAX_FACES * 3;
    }
    if (remainder != 0)
        glDrawArrays(GL_TRIANGLES, start_index, remainder * 3);

    rendering_cleanup();
}

}