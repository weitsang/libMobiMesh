/*
 *  GLVpmMeshRenderer.cpp
 *  MobiMesh
 *
 *  Created by Ong Yuh Shin on 3/29/11.
 *  Copyright 2011 National University of Singapore. All rights reserved.
 *
 */

#include "GLVpmMeshRenderer.h"
#include "ViewingParameters.h"

namespace MobiMesh {

GLVpmMeshRenderer::GLVpmMeshRenderer(GLMeshRenderer& renderer,
                                     ViewDependentProgressiveMesh* mesh)
	: GLViewDependentMeshRenderer(renderer), mesh_(mesh),
    need_adapt_(true)
{
    set_mesh(mesh_);
    set_default_view();
    activate();
}

GLVpmMeshRenderer::~GLVpmMeshRenderer() { }

void GLVpmMeshRenderer::activate()
{
    if (mesh_)
        mesh_->activate();
    renderer_.activate();
}

void GLVpmMeshRenderer::deactivate()
{
    renderer_.deactivate();
    if (mesh_)
        mesh_->deactivate();
}

bool GLVpmMeshRenderer::adapt_mesh()
{
    if (!mesh_)
        return false; // no mesh to render

    if (view_changed())
    {
        // view has changed, read in the new modelview matrix
        float* modelview_matrix = viewing_params_.modelview_matrix();

        glm::mat4 lookat_matrix, transformation_matrix;
        capture_view(lookat_matrix, transformation_matrix);
        glm::mat4 current_modelview_matrix =
            lookat_matrix * transformation_matrix;

        for (int col = 0; col < 4; ++col)
        {
            for (int row = 0; row < 4; ++row)
            {
                modelview_matrix[col*4 + row] =
                    current_modelview_matrix[col][row];
            }
        }

        viewing_params_.set_aspect((float) get_viewport_width() /
                                   (float) get_viewport_height());
        viewing_params_.update_viewing_configurations();

        // Refine fewer vertices to improve framerate (as user is changing
        // views) while still improving quality.
        mesh_->set_viewing_parameters(viewing_params_);
    }

    mesh_->adapt(MAX_REFINEMENTS);

    return true;
}

}
