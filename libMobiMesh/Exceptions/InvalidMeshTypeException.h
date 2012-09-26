/*
 *  InvalidMeshTypeException.h
 *  MobiMesh
 *
 *  Created by Ong Yuh Shin on 3/28/11.
 *  Copyright 2011 National University of Singapore. All rights reserved.
 *
 */

#ifndef __INVALID_MESH_TYPE_EXCEPTION_H_
#define __INVALID_MESH_TYPE_EXCEPTION_H_

#include "MobiMeshException.h"

namespace MobiMesh {

class InvalidMeshTypeException : public MobiMeshException
{
public:
	InvalidMeshTypeException(const std::string& err_msg="")
		: MobiMeshException(err_msg) { }

	virtual ~InvalidMeshTypeException() throw() { }

    const char* get_exception_name() const throw() { return "InvalidMeshTypeException"; }
};
	
}

#endif