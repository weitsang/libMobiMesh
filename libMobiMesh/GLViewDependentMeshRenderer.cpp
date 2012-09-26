//
//  GLViewDependentMeshRenderer.cpp
//  MobiMesh
//
//  Created by Ong Yuh Shin on 5/19/11.
//  Copyright 2011 National University of Singapore. All rights reserved.
//

#include "GLViewDependentMeshRenderer.h"
#include "InactiveRendererException.h"

namespace MobiMesh {

GLViewDependentMeshRenderer::
GLViewDependentMeshRenderer(GLMeshRenderer& renderer)
    : renderer_(renderer), view_changed_(true) { }

void GLViewDependentMeshRenderer::render()
{
    // Only update the vbos when the when the view has not changed
    // since the previous frame.
    if (adapt_mesh())
        renderer_.update_mesh();

    renderer_.render();
    view_changed_ = false;
}

}