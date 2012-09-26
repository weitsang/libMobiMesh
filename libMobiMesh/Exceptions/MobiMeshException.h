/*
 *  MobiMeshException.h
 *  MobiMesh
 *
 *  Created by Ong Yuh Shin on 3/14/11.
 *  Copyright 2011 National University of Singapore. All rights reserved.
 *
 */

#ifndef __MOBI_MESH_EXCEPTION_H_
#define __MOBI_MESH_EXCEPTION_H_

#include <exception>
#include <string>

namespace MobiMesh {

/**
 * Base exception class for all exceptions thrown by the MobiMesh API.
 */
class MobiMeshException : public std::exception
{
public:
	MobiMeshException(const std::string& err_msg="") : err_msg_(err_msg) { }
	virtual ~MobiMeshException() throw() { }

	virtual const char* what() const throw()
    {
        std::string what_str(get_exception_name());
        what_str += ": ";
        what_str += err_msg_;
        return what_str.c_str();
    }

    virtual const char* get_exception_name() const throw() { return "MobiMeshException"; }

private:
	std::string err_msg_;
};

}

#endif