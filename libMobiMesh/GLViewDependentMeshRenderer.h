//
//  GLViewDependentMeshRenderer.h
//  MobiMesh
//
//  Created by Ong Yuh Shin on 5/19/11.
//  Copyright 2011 National University of Singapore. All rights reserved.
//

#ifndef __GL_VIEW_DEPENDENT_MESH_RENDERER_H_
#define __GL_VIEW_DEPENDENT_MESH_RENDERER_H_

#include "GLMeshRenderer.h"

namespace MobiMesh {

/**
 * Abstract base class for rendering one of the view-dependent meshes (.vpm,
 * .tvpm, .tvpm64)
 * Concrete classes will implement the algorithm required for adapting the
 * respective mesh being rendered. This algorithm is called whenever the
 * render() method of this class is called.
 */
class GLViewDependentMeshRenderer : public IMeshRenderer
{
public:

    ///////////// Constructors

    /**
     * Creates a renderer for rendering one of the view dependent meshes.
     *
     * @param[in] renderer The renderer to be used for doing the actual mesh
     * rendering using the appropriate version of OpenGL ES. Caller must
     * ensure that lifetime of supplied renderer outlives that of
     * the created renderer.
     */
    explicit GLViewDependentMeshRenderer(GLMeshRenderer& renderer);

    virtual ~GLViewDependentMeshRenderer() { }

    ///////////// Methods inherited from GLMeshRenderer

    void render();

    void translate(float dx, float dy, float dz)
    {
        view_changed_ = true;
        renderer_.translate(dx, dy, dz);
    }

    void rotate(float angle_x_radians, float axis_x, float axis_y, float axis_z)
    {
        view_changed_ = true;
        renderer_.rotate(angle_x_radians, axis_x, axis_y, axis_z);
    }

    void scale(float scale_factor)
    {
        view_changed_ = true;
        renderer_.scale(scale_factor);
    }

    void set_default_view()
    {
        view_changed_ = true;
        renderer_.set_default_view();
    }

    void capture_view(glm::mat4& lookat, glm::mat4& transformation) const
    { renderer_.capture_view(lookat, transformation); }

    void set_view(const glm::mat4& lookat, const glm::mat4 transformation)
    {
        view_changed_ = true;
        renderer_.set_view(lookat, transformation);
    }

    float get_bounding_length() const
    { return renderer_.get_bounding_length(); }

    int get_viewport_width() const { return renderer_.get_viewport_width(); }
    int get_viewport_height() const { return renderer_.get_viewport_height(); }

    virtual void activate() =0;
    virtual void deactivate() =0;
    bool is_active() { return renderer_.is_active(); }

    /// Max number of refinements per frame
    static const int MAX_REFINEMENTS = 1000;

protected:

    /**
     * Refines mesh based on current view. Subclasses must implement this
     * method to determine how the mesh is to be rendered.
     * @return true iff the mesh was modified/adapted
     */
    virtual bool adapt_mesh() =0;

    /**
     * For subclasses to set the mesh to be rendered by the renderer.
     */
    void set_mesh(IMesh* mesh) { renderer_.set_mesh(mesh); }

    bool view_changed() { return view_changed_; }

    GLMeshRenderer& renderer_;

private:
    
    GLViewDependentMeshRenderer(const GLViewDependentMeshRenderer&);
    GLViewDependentMeshRenderer& operator=(const GLViewDependentMeshRenderer&);
    
    bool view_changed_;
};

}

#endif