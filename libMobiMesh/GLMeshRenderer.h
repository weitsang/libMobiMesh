//
//  GLMeshRenderer.h
//  MobiMesh
//
//  Created by Ong Yuh Shin on 5/5/11.
//  Copyright 2011 National University of Singapore. All rights reserved.
//

#ifndef __GL_MESH_RENDERER_H_
#define __GL_MESH_RENDERER_H_

#include <vector>

#include <OpenGLES/ES1/gl.h>
#include <OpenGLES/ES1/glext.h>
#include <OpenGLES/ES2/gl.h>

// #include "matrix_transform.hpp"

#include "GLMeshRendererViewingParameters.h"
#include "IMesh.h"
#include "IMeshRenderer.h"

namespace MobiMesh {

/**
 * Abstract base class for mesh rendering using OpenGL.
 * Subclasses of this must override the pure virtual functions in this class
 * using specific OpenGL versions/rendering pipelines.
 *
 * The render method cannot be overridden. The rendering process is
 * fixed. The GLMeshRenderer uses VBOs containing the vertices, normals and
 * faces respectively, and is capable of rendering any mesh that implements
 * the IMesh interface.
 */
class GLMeshRenderer : public IMeshRenderer
{
public:
    
	///////////// Constructors

    /**
     * Create a renderer with a viewport size defined by width and height.
     * @param[in] width viewport width
     * @param[in] height viewport height
     */
	GLMeshRenderer(GLsizei width, GLsizei height);

    /**
     * Create a renderer with a viewport size defined by width and height.
     * @param[in] mesh caller supplied mesh. Caller needs to ensure that the
     * lifetime of the supplied object exceeds that of the created
     * GLMeshRenderer object.
     * @param[in] width viewport width
     * @param[in] height viewport height
     */
	GLMeshRenderer(IMesh* mesh, GLsizei width, GLsizei height);

	virtual ~GLMeshRenderer();
    
	///////////// Methods inherited from IMeshRenderer

    /**
     * @throw InvalidGeometryException if the mesh geometry is incorrect
     * (and therefore unable to render)
     */
    virtual void render();

	virtual void translate(float dx, float dy, float dz)
    { view_params_.translate(dx, dy, dz); }
    
	virtual void rotate(float angle_degrees,
                        float axis_x, float axis_y, float axis_z)
    { view_params_.rotate(angle_degrees, axis_x, axis_y, axis_z); }
    
    virtual void scale(float scale_factor)
    { view_params_.scale(scale_factor); }

    void set_default_view() { view_params_.set_default_view(); }

    void capture_view(glm::mat4& lookat, glm::mat4& transformation) const
    {
        lookat = view_params_.get_lookat_matrix();
        transformation = view_params_.get_transformation_matrix();
    }

    void set_view(const glm::mat4& lookat, const glm::mat4 transformation)
    {
        view_params_.set_lookat_matrix(lookat);
        view_params_.set_transformation_matrix(transformation);
    }

    float get_bounding_length() const
    { return view_params_.get_bounding_length(); }

    int get_viewport_width() const { return viewport_width_; }
    int get_viewport_height() const { return viewport_height_; }

    /**
     * Activates the mesh for rendering.
     * Note that the caller needs to re-supply the color buffer using
     * set_colors if the mesh was previously de-activated and if the
     * caller wants to render with colors.
     */
    void activate();

    /**
     * Frees up OpenGL buffers. Mesh cannot be rendered in this state.
     */
    void deactivate();

    bool is_active() { return is_active_; }

    ///////////// Methods access/update rendering params

    /**
     * This should be called when the underlying mesh to be rendered has
     * changed, i.e. the points/normals/faces have changed. This is so that
     * the renderer can update its vbo.
     */
    void update_mesh();

    void set_mesh(IMesh* mesh)
    { view_params_.set_mesh(mesh); update_mesh(); }

    /**
     * @return true iff the renderer will but supplying a color pointer
     * to OpenGL when rendering.
     * set_color must be called to set the color array if this function
     * returns true.
     */
    bool render_with_color() const { return with_color_; }

    /**
     * @param[in] true iff a color pointer should be used when rendering.
     * If rendering with color, the color vector should be supplied using
     * set_color.
     */
    void render_with_color(bool color); 
    /**
     * Sets the color array to be used when rendering.
     */
    void set_color(const std::vector<OpenMesh::Vec4uc>& color);

protected:

    // Recommended FOVY angle in degrees to use
    static const int FOVY = 45;

    GLMeshRendererViewingParameters view_params_;

    ///////////// Abstract methods: Subclasses must implement these because
    ///////////// the specific ways to perform such setups are OpenGL
    ///////////// version specific.
    
    /**
     * Sets up the modelview and projection matrices accordingly
     * i.e. lookat, translate, scale, rotate, etc.
     */
    virtual void setup_modelview_projection() =0;

    /**
     * Called to tell OpenGL which vertices need to be rendered.
     * Prior to calling this method, the corresponding buffer containing all
     * mesh vertices will have been bounded.
     * The vertices are 3D coordinates of type float
     *
     * @param[in] pointer the byte offset into the buffer object's data
     * store, referring to the first vertex in the array where we want to
     * draw from.
     */
    virtual void setup_vertices(const GLvoid* pointer) =0;

    /**
     * Called to tell OpenGL which normals need to be used.
     * Prior to calling this method, the corresponding buffer containing all
     * vertex normals will have been bounded.
     * The vertices are 3D coordinates of type float
     *
     * @param[in] pointer the byte offset into the buffer object's data
     * store, referring to the first vertex normal in the array where we
     * want to draw from.
     */
    virtual void setup_normals(const GLvoid* pointer) =0;

    /**
     * Called to tell OpenGL which colors need to be used
     * Prior to calling this method, the corresponding buffer containing all
     * colors will have been bounded.
     * The vertices are 3D coordinates of type float
     *
     * @param[in] pointer byte offset into the buffer object's data
     * store, referring tot he first color element in the array where
     * we want to draw from.
     */
    virtual void setup_colors(const GLvoid* pointer) =0;

    /**
     * Called before rendering is done. Calls to enable client state or
     * vertex attribute arrays should call here.
     */
    virtual void rendering_setup() =0;

    /**
     * Called once rendering is done. Calls to disable client state or vertex attribute arrays
     * should be implemented here.
     */
    virtual void rendering_cleanup() =0;

private:

    GLMeshRenderer(const GLMeshRenderer&);
    GLMeshRenderer& operator=(const GLMeshRenderer&);

    GLuint vbo_vertices_, vbo_normals_, vbo_faces_, vbo_colors_;

    // Flag to determine rendering with color (true) or with normals (false)
    bool with_color_;

    bool is_active_;

    // Viewport size
    GLsizei viewport_width_;
    GLsizei viewport_height_;

    // Number of faces in VBO
    size_t num_faces_;

    // view_params_.get_mesh()->has_faces() must be true
    void bind_and_draw_vbo_draw_arrays();

    // view_params_.get_mesh()->has_face_vertices() must be true
    void bind_and_draw_vbo_draw_elements();
};

}

#endif