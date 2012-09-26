//
//  InvalidGeometryException.h
//  MobiMesh
//
//  Created by Ong Yuh Shin on 5/26/11.
//  Copyright 2011 National University of Singapore. All rights reserved.
//

#ifndef __INVALID_GEOMETRY_EXCEPTION_H_
#define __INVALID_GEOMETRY_EXCEPTION_H_

#include "MobiMeshException.h"

namespace MobiMesh {
    
/**
 * Exception class indicating an error in the geometry of the mesh
 */
class InvalidGeometryException : public MobiMeshException
{
public:

    InvalidGeometryException(const std::string& err_msg="") 
        : MobiMeshException(err_msg) { }

    virtual ~InvalidGeometryException() throw() { }
    
    const char* get_exception_name() const throw() { return "InvalidGeometryException"; }
};
    
}

#endif