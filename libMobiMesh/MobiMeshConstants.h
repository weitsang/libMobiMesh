//
//  MobiMeshConstants.h
//  MobiMesh
//
//  Created by Ong Yuh Shin on 6/9/11.
//  Copyright 2011 National University of Singapore. All rights reserved.
//

#ifndef __MOBI_MESH_CONSTANTS_H_
#define __MOBI_MESH_CONSTANTS_H_

namespace MobiMesh {
namespace Constants {

/// Max number of unique vertices supported by the renderers that use the
/// glDrawElements call. This is due to OpenGL ES's limitation on the size
/// of the vertex indices (GL_UNSIGNED_SHORT)
const int MAX_VERTICES = 0xFFFF;

/// Max number of faces allowed when using glDrawArrays with triangle
/// primitives. This needs to be limited because glDrawArrays fails if there
/// are too many faces to draw.
const unsigned int MAX_FACES = 0xFFFF;

}
}

#endif