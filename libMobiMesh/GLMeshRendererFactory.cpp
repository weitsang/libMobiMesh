/*
 *  GLMeshRendererFactory.cpp
 *  MobiMesh
 *
 *  Created by Ong Yuh Shin on 3/28/11.
 *  Copyright 2011 National University of Singapore. All rights reserved.
 *
 */

#include "GLMeshRendererFactory.h"
#include "ProgressiveMesh.h"
#include "ViewDependentProgressiveMesh.h"
#include "TrulyViewDependentPM32.h"
#include "TrulyViewDependentPM64.h"
#include "GLMeshRendererES1.h"
#include "GLVpmMeshRenderer.h"
#include "GLTvpmMeshRenderer.h"
#include "GLMeshRendererES2.h"

#include "InvalidMeshTypeException.h"
#include "UnableToCreateRendererException.h"

namespace MobiMesh {

GLMeshRendererFactory::GLMeshRendererFactory()
    : renderers_(MAX_RENDERERS),
    renderers_for_vpm_(MAX_RENDERERS),
    meshes_(MAX_RENDERERS)
{
    for (int i=0; i<MAX_RENDERERS; ++i)
    {
        renderers_[i] = NULL;
        renderers_for_vpm_[i] = NULL;
        meshes_[i] = NULL;
        free_handles_.push_back(i);
    }
}

GLMeshRendererFactory::~GLMeshRendererFactory()
{
    for (std::vector<IMeshRenderer*>::iterator renderer = renderers_.begin();
         renderer != renderers_.end(); ++renderer)
    {
        if (*renderer)
            delete *renderer;
    }

    for (std::vector<GLMeshRenderer*>::iterator
         renderer = renderers_for_vpm_.begin();
         renderer != renderers_for_vpm_.end();
         ++renderer)
    {
        if (*renderer)
            delete *renderer;
    }

    for (std::vector<IMesh*>::iterator mesh = meshes_.begin();
         mesh != meshes_.end(); ++mesh)
    {
        if (*mesh)
            delete *mesh;
    }
}

GLMeshRendererFactory::RendererHandle
GLMeshRendererFactory::new_mesh_renderer(const std::string& filename,
                                         GLsizei viewport_width,
                                         GLsizei viewport_height)
{
    // Get an available handle
    RendererHandle handle = get_next_avail_handle();

    // No more free handles, do not create the object
    if (!is_valid_handle(handle))
        return handle;

	if (filename.rfind(".vpm") != std::string::npos)
	{
        ViewDependentProgressiveMesh* mesh = NULL;
        GLMeshRendererES1* renderer = NULL;
        GLVpmMeshRenderer* result = NULL;
        try
        {
            mesh = new ViewDependentProgressiveMesh(filename);
            renderer = new GLMeshRendererES1(viewport_width, viewport_height);
            result = new GLVpmMeshRenderer(*renderer, mesh);
        }
        catch (MobiMeshException& ex)
        {
            if (mesh)
                delete mesh;
            if (renderer)
                delete renderer;
            if (result)
                delete result;
            throw UnableToCreateRendererException(ex.what());
        }
        renderers_[handle] = result;
        renderers_for_vpm_[handle] = renderer;
        meshes_[handle] = mesh;
	}
	else if (filename.rfind( ".pm") != std::string::npos)
	{
        ProgressiveMesh* mesh = NULL;
        GLMeshRendererES1* renderer = NULL;
        try
        {
            mesh = new ProgressiveMesh(filename);
            renderer = new GLMeshRendererES1(mesh,
                                             viewport_width, viewport_height);
        }
        catch (MobiMeshException& ex)
        {
            if (mesh)
                delete mesh;
            if (renderer)
                delete renderer;
            throw UnableToCreateRendererException(ex.what());
        }
        renderers_[handle] = renderer;
        meshes_[handle] = mesh;
	}
    else if (filename.rfind(".tvpm") != std::string::npos)
    {
        TrulyViewDependentPM* mesh = NULL;
        GLMeshRendererES1* renderer = NULL;
        GLTvpmMeshRenderer* result = NULL;
        try
        {
            if (filename.rfind(".tvpm64") != std::string::npos)
                mesh = new TrulyViewDependentPM64(filename);
            else
                mesh = new TrulyViewDependentPM32(filename);
            renderer = new GLMeshRendererES1(viewport_width, viewport_height);
            result = new GLTvpmMeshRenderer(*renderer, mesh);
        }
        catch (MobiMeshException& ex)
        {
            if (mesh)
                delete mesh;
            if (renderer)
                delete renderer;
            if (result)
                delete result;
            throw UnableToCreateRendererException(ex.what());
        }
        renderers_[handle] = result;
        renderers_for_vpm_[handle] = renderer;
        meshes_[handle] = mesh;
    }
	else
    {
		throw InvalidMeshTypeException("filename must end with .vpm, .pm, "
                                       ".tvpm, or .tvpm64");
    }

    use_next_avail_handle();
	return handle;
}

GLMeshRendererFactory::RendererHandle
GLMeshRendererFactory::new_es2_mesh_renderer(
    const std::string& mesh_filename,
    const char* vertex_shader,
    const char* fragment_shader,
    GLsizei viewport_width,
    GLsizei viewport_height)
{
    // Get an available handle
    RendererHandle handle = get_next_avail_handle();

    // No more free handles, do not create the object
    if (!is_valid_handle(handle))
        return handle;

    if (mesh_filename.rfind(".vpm") != std::string::npos)
    {
        ViewDependentProgressiveMesh* mesh = NULL;
        GLMeshRendererES2* renderer = NULL;
        GLVpmMeshRenderer* result = NULL;
        try
        {
            mesh = new ViewDependentProgressiveMesh(mesh_filename);
            renderer = new GLMeshRendererES2(vertex_shader, fragment_shader,
                                             viewport_width, viewport_height);
            result = new GLVpmMeshRenderer(*renderer, mesh);
        }
        catch (MobiMeshException& ex)
        {
            if (mesh)
                delete mesh;
            if (renderer)
                delete renderer;
            if (result)
                delete result;
            throw UnableToCreateRendererException(ex.what());
        }
        renderers_[handle] = result;
        renderers_for_vpm_[handle] = renderer;
        meshes_[handle] = mesh;
    }
    else if (mesh_filename.rfind( ".pm") != std::string::npos)
    {
        ProgressiveMesh* mesh = NULL;
        GLMeshRendererES2* renderer = NULL;
        try
        {
            mesh = new ProgressiveMesh(mesh_filename);
            renderer = new GLMeshRendererES2(vertex_shader, fragment_shader,
                                             mesh,
                                             viewport_width, viewport_height);
        }
        catch (MobiMeshException& ex)
        {
            if (mesh)
                delete mesh;
            if (renderer)
                delete renderer;
            throw UnableToCreateRendererException(ex.what());
        }
        renderers_[handle] = renderer;
        meshes_[handle] = mesh;
    }
    else if (mesh_filename.rfind(".tvpm") != std::string::npos)
    {
        TrulyViewDependentPM* mesh = NULL;
        GLMeshRendererES2* renderer = NULL;
        GLTvpmMeshRenderer* result = NULL;
        try
        {
            if (mesh_filename.rfind(".tvpm64") != std::string::npos)
                mesh = new TrulyViewDependentPM64(mesh_filename);
            else
                mesh = new TrulyViewDependentPM32(mesh_filename);
            renderer = new GLMeshRendererES2(vertex_shader, fragment_shader,
                                             mesh,
                                             viewport_width, viewport_height);
            result = new GLTvpmMeshRenderer(*renderer, mesh);
        }
        catch (MobiMeshException& ex)
        {
            if (mesh)
                delete mesh;
            if (renderer)
                delete renderer;
            if (result)
                delete result;
            throw UnableToCreateRendererException(ex.what());
        }
        renderers_[handle] = result;
        renderers_for_vpm_[handle] = renderer;
        meshes_[handle] = mesh;
    }
	else
    {
		throw InvalidMeshTypeException("filename must end with .vpm, .pm, "
                                       ".tvpm, or .tvpm64");
    }

    // If everything was successful, i.e. no exception thrown, mark this
    // handle as used
    use_next_avail_handle();

    return handle;
}

void GLMeshRendererFactory::free_mesh_renderer(RendererHandle handle)
{
    // Invalid handle, nothing to free
    if (!is_valid_handle(handle))
        return;

    if (renderers_[handle]) // valid non-null renderer, needs to be freed
    {
        delete renderers_[handle];
        renderers_[handle] = NULL;

        if (renderers_for_vpm_[handle])
        {
            delete renderers_for_vpm_[handle];
            renderers_for_vpm_[handle] = NULL;
        }

        if (meshes_[handle])
        {
            delete meshes_[handle];
            meshes_[handle] = NULL;
        }

        free_handles_.push_back(handle);
    }
}

IMeshRenderer* GLMeshRendererFactory::get_renderer(RendererHandle handle)
{
    if (!is_valid_handle(handle))
        return NULL;
    else
        return renderers_[handle];
}

}
