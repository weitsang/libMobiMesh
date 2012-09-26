/*
 *  VHierarchyNodeIndex.h
 *  MobiMesh
 *
 *  Created by Ong Yuh Shin on 3/10/11.
 *  Copyright 2011 National University of Singapore. All rights reserved.
 *
 */

#ifndef __VHIERARCHY_NODE_INDEX_H_
#define __VHIERARCHY_NODE_INDEX_H_

#include "VectorT.h"

namespace MobiMesh {

/**
 * Represents a node index in a view dependent progressive mesh
 */
class VHierarchyNodeIndex
{

public:
	//////////////// Constructors
    VHierarchyNodeIndex() : value_(0) { }

	explicit VHierarchyNodeIndex(unsigned int _value)
		: value_(_value) { }

	VHierarchyNodeIndex(const VHierarchyNodeIndex &_other)
		: value_(_other.value_) { }

    VHierarchyNodeIndex(
	  unsigned int _tree_id, unsigned int _node_id, unsigned short _tree_id_bits)
	{
		assert(_tree_id < ((unsigned int) 0x00000001 << _tree_id_bits));
		assert(_node_id < ((unsigned int) 0x00000001 << (32 - _tree_id_bits)));
		value_ = (_tree_id << (32 - _tree_id_bits)) | _node_id;
	}

	///////////// Misc

	bool is_valid(unsigned short _tree_id_bits) const
	{
		return node_id(_tree_id_bits) != 0 ? true : false;
	}

	unsigned int tree_id(unsigned short _tree_id_bits) const
	{
		return value_ >> (32 - _tree_id_bits);
	}

    unsigned int node_id(unsigned short _tree_id_bits) const
	{
		return value_ & ((unsigned int) 0xFFFFFFFF >> _tree_id_bits);
	}

	bool operator< (const VHierarchyNodeIndex &other) const
	{
		return (value_ < other.value_) ? true : false;
	}

	unsigned int value() const
	{
		return value_;
	}

private:
	unsigned int value_;

};

}

#endif