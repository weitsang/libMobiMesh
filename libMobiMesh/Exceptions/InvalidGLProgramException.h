//
//  InvalidGLProgramException.h
//  MobiMesh
//
//  Created by Ong Yuh Shin on 4/11/11.
//  Copyright 2011 National University of Singapore. All rights reserved.
//

#ifndef __INVALID_GL_PROGRAM_EXCEPTION_H_
#define __INVALID_GL_PROGRAM_EXCEPTION_H_

#include "MobiMeshException.h"

namespace MobiMesh {

/**
 * Exception class indicating an error validating a program in
 * OpenGL ES 2.0.
 */
class InvalidGLProgramException : public MobiMeshException
{
public:

    InvalidGLProgramException(const std::string& err_msg="") 
        : MobiMeshException(err_msg) { }

    virtual ~InvalidGLProgramException() throw() { }
    
    const char* get_exception_name() const throw() { return "InvalidGLProgramException"; }
};

}

#endif