/*
 *  InvalidFileFormat.h
 *  MobiMesh
 *
 *  Created by Ong Yuh Shin on 3/14/11.
 *  Copyright 2011 National University of Singapore. All rights reserved.
 *
 */

#ifndef __INVALID_FILE_FORMAT_EXCEPTION_H_
#define __INVALID_FILE_FORMAT_EXCEPTION_H_

#include <string>
#include "MobiMeshException.h"

namespace MobiMesh {

class InvalidFileFormatException : public MobiMeshException
{
public:
	InvalidFileFormatException(const std::string& err_msg="")
	  : MobiMeshException(err_msg) { }

	/**
	 * @param[in] expected_format the file format that is expected
	 * @param[in] supplied_format the file format that was received
	 */
	InvalidFileFormatException(
	  const std::string& err_msg,
	  const std::string& expected_format,
	  const std::string& supplied_format) :
		MobiMeshException(err_msg),
		expected_format_(expected_format),
		supplied_format_(supplied_format) { }

	virtual ~InvalidFileFormatException() throw() { }

	std::string get_expected_format() { return expected_format_; }
	std::string get_supplied_format() { return supplied_format_; }

private:
	std::string expected_format_;
	std::string supplied_format_;
};

}

#endif
