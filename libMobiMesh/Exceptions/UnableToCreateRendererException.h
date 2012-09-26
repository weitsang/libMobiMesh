//
//  UnableToCreateRendererException.h
//  MobiMesh
//
//  Created by Ong Yuh Shin on 5/24/11.
//  Copyright 2011 National University of Singapore. All rights reserved.
//

#ifndef __UNABLE_TO_CREATE_RENDERER_EXCEPTION_H_
#define __UNABLE_TO_CREATE_RENDERER_EXCEPTION_H_

#include "MobiMeshException.h"

namespace MobiMesh {
    
/**
 * Exception class indicating an error creating a renderer
 */
class UnableToCreateRendererException : public MobiMeshException
{
public:

    UnableToCreateRendererException(const std::string& err_msg="") 
        : MobiMeshException(err_msg) { }
        
    virtual ~UnableToCreateRendererException() throw() { }
    
    const char* get_exception_name() const throw() { return "UnableToCreateRendererException"; }
};

}

#endif