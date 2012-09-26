//
//  ProgramLinkingException.h
//  MobiMesh
//
//  Created by Ong Yuh Shin on 4/10/11.
//  Copyright 2011 National University of Singapore. All rights reserved.
//

#ifndef __PROGRAM_LINKING_EXCEPTION_H_
#define __PROGRAM_LINKING_EXCEPTION_H_

#include "MobiMeshException.h"

namespace MobiMesh {
    
/**
 * Exception class indicating an error linking a program in
 * OpenGL ES 2.0.
 */
class ProgramLinkingException : public MobiMeshException
{
public:

    ProgramLinkingException(const std::string& err_msg="") 
        : MobiMeshException(err_msg) { }

    virtual ~ProgramLinkingException() throw() { }
    
    const char* get_exception_name() const throw() { return "ProgramLinkingException"; }
};
    
}

#endif