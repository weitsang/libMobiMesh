/*
 *  GLMeshRendererFactory.h
 *  MobiMesh
 *
 *  Created by Ong Yuh Shin on 3/28/11.
 *  Copyright 2011 National University of Singapore. All rights reserved.
 *
 */

#ifndef __GL_MESH_RENDERER_FACTORY_H_
#define __GL_MESH_RENDERER_FACTORY_H_

// for GLsizei type
#include <OpenGLES/ES2/gl.h>

#include <list>
#include <vector>

#include "IMeshRenderer.h"
#include "GLMeshRenderer.h"
#include "GLVpmMeshRenderer.h"

namespace MobiMesh {

/**
 * A factory class responsible for creating GLMeshRenderers based on the
 * supplied file type. It is also responsible for managing memory of the
 * created renderer objects.
 */
class GLMeshRendererFactory {

public:

    GLMeshRendererFactory();
    ~GLMeshRendererFactory();

    /// Maximum number of active renders at any point in time.
    static const int MAX_RENDERERS = 100;

    typedef int RendererHandle;
    static const RendererHandle INVALID_HANDLE = -1;

    /**
     * Loads the mesh from the supplied filename. Mesh type is determined by
     * filename extension. The appropriate renderer for the mesh type is
     * returned. The returned handle is guaranteed to be a previously unused
     * handle at the time where this method is called. However, it might
     * reuse previously freed handles.
     *
     * @param[in] mesh_filename filename of file containing mesh data
     * @return RendererHandle for the mesh renderer, or an invalid handle
     * if the factory cannot support any more renderers.
     * @throw InvalidMeshTypeException if filename does not end with .pm
     * or .vpm
     * @throw UnableToCreateRendererException if there was an error creating
     * the renderer
     */
    RendererHandle new_mesh_renderer(const std::string& mesh_filename,
                                     GLsizei viewport_width,
                                     GLsizei viewport_height);

    /**
     * Same as get_mesh_renderer(const std::string&) but for renderers using
     * OpenGL ES2.0
     * @throw InvalidMeshTypeException if filename does not end with .pm
     * or .vpm
     * @throw UnableToCreateRendererException if there was an error creating
     * the renderer
     */
    RendererHandle new_es2_mesh_renderer(const std::string& mesh_filename,
                                         const char* vertex_shader,
                                         const char* fragment_shader,
                                         GLsizei viewport_width,
                                         GLsizei viewport_height);

    /**
     * Free up memory allocated to the renderer. No effect if handle is
     * invalid.
     */
    void free_mesh_renderer(RendererHandle handle);

    /**
     * @return pointer to the renderer with the given handle, or NULL if
     * the handle does not refer to a valid handle
     */
    IMeshRenderer* get_renderer(RendererHandle handle);

private:

    GLMeshRendererFactory(const GLMeshRendererFactory&);
    GLMeshRendererFactory& operator=(const GLMeshRendererFactory&);

    // List of renderers created and returned to the caller
    // The renderer's handle is the index to the element that contains it.
    // A NULL element indicates a free slot.
    std::vector<IMeshRenderer*> renderers_;

    // List of renderers created to support GLVpmMeshRenderer creation,
    // which requires that a reference to a GLMeshRenderer be passed in.
    // The GLVpmMeshRenderer object's handle is the index to its
    // corresponding GLMeshRenderer object in the vector.
    // A NULL element indicates a free slot and that the index/handle at
    // that slot is not a handle of a vpm renderer.
    std::vector<GLMeshRenderer*> renderers_for_vpm_;

    // List of meshes for each renderer
    std::vector<IMesh*> meshes_;

    // List of unused handles.
    std::list<RendererHandle> free_handles_;

    bool is_valid_handle(RendererHandle handle)
    { return handle >= 0 && handle < MAX_RENDERERS; }

    // Gets the next available handle (handle is still avail after this call)
    RendererHandle get_next_avail_handle() const
    {
        if (free_handles_.size() == 0)
            return -1;

        RendererHandle free_handle = free_handles_.front();
        return free_handle;
    }

    // Uses the next available handle (handle no longer avail after this call)
    RendererHandle use_next_avail_handle()
    {
        if (free_handles_.size() == 0)
            return -1;
        
        RendererHandle free_handle = free_handles_.front();
        free_handles_.pop_front();
        return free_handle;
    }
};

}

#endif
