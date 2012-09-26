//
//  InvalidThresholdException.h
//  MobiMesh
//
//  Created by Ong Yuh Shin on 6/9/11.
//  Copyright 2011 National University of Singapore. All rights reserved.
//

#ifndef __INVALID_THRESHOLD_EXCEPTION_H_
#define __INVALID_THRESHOLD_EXCEPTION_H_

#include "MobiMeshException.h"

namespace MobiMesh {
    
/**
 * Exception class indicating invalid threshold values for TVPM mesh
 * refinement
 */
class InvalidThresholdException : public MobiMeshException
{
public:

    InvalidThresholdException(const std::string& err_msg="") 
        : MobiMeshException(err_msg) { }

    virtual ~InvalidThresholdException() throw() { }

    const char* get_exception_name() const throw()
    { return "InvalidThresholdException"; }
};
    
}

#endif