//
//  InactiveRendererException.h
//  MobiMesh
//
//  Created by Ong Yuh Shin on 6/10/11.
//  Copyright 2011 National University of Singapore. All rights reserved.
//

#ifndef __INACTIVE_RENDERER_EXCEPTION_H_
#define __INACTIVE_RENDERER_EXCEPTION_H_

#include "MobiMeshException.h"

namespace MobiMesh {
    
/**
 * Exception class indicating that a method other than activate() or
 * deactivate() was called on a renderer object that is in inactive state.
 */
class InactiveRendererException : public MobiMeshException
{
public:

    InactiveRendererException(const std::string& err_msg="") 
        : MobiMeshException(err_msg) { }

    virtual ~InactiveRendererException() throw() { }

    const char* get_exception_name() const throw()
    { return "InactiveRendererException"; }

};
    
}

#endif