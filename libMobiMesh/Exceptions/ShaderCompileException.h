//
//  ShaderCompileException.h
//  MobiMesh
//
//  Created by Ong Yuh Shin on 4/10/11.
//  Copyright 2011 National University of Singapore. All rights reserved.
//

#ifndef __SHADER_COMPILE_EXCEPTION_H_
#define __SHADER_COMPILE_EXCEPTION_H_

#include "MobiMeshException.h"

namespace MobiMesh {
    
/**
 * Exception class indicating an error compiling a shader in
 * OpenGL ES 2.0.
 */
class ShaderCompileException : public MobiMeshException
{
public:

    ShaderCompileException(
        const std::string& shader_src, const std::string& err_msg="") 
        : MobiMeshException(err_msg), shader_src_(shader_src) { }

    virtual ~ShaderCompileException() throw() { }
    
    std::string get_shader_src() { return shader_src_; }
    
    const char* get_exception_name() const throw() { return "ShaderCompileException"; }

private:
    std::string shader_src_;
};

}

#endif