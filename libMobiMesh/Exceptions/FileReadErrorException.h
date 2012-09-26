/*
 *  FileReadErrorException.h
 *  MobiMesh
 *
 *  Created by Ong Yuh Shin on 3/24/11.
 *  Copyright 2011 National University of Singapore. All rights reserved.
 *
 */

#ifndef __FILE_READ_ERROR_EXCEPTION_H_
#define __FILE_READ_ERROR_EXCEPTION_H_

#include "MobiMeshException.h"

namespace MobiMesh {

class FileReadErrorException : public MobiMeshException
{
public:
	FileReadErrorException(const std::string& err_msg="")
		: MobiMeshException(err_msg) { }

	virtual ~FileReadErrorException() throw() { }

    const char* get_exception_name() const throw() { return "FileReadErrorException"; }
};
	
}

#endif
