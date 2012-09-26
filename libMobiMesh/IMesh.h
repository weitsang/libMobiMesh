/*
 *  IMesh.h
 *  MobiMesh
 *
 *  Created by Ong Yuh Shin on 3/8/11.
 *  Copyright 2011 National University of Singapore. All rights reserved.
 *
 */

#ifndef __IMESH_H_
#define __IMESH_H_

#include <string>
#include <vector>

#include "VectorT.h"

namespace MobiMesh {

class IMeshRenderer;

/**
 * Abstract base class/interface for all meshes in the library.
 */
class IMesh
{
public:

	///////////// Types and constants

	typedef OpenMesh::Vec3f Point;
	typedef OpenMesh::Vec3f Normal;
	typedef std::vector<Point> PointContainer;
	typedef std::vector<Normal> VertexNormalContainer;
	typedef std::vector<unsigned short> FaceContainer;

	///////////// Constructors

	virtual ~IMesh() { }
	
	///////////// Loading mesh

	/**
	 * Loads a mesh from file, overwriting mesh data that is currently
	 * stored in this object (if any)
	 * @param[in] filename name of file containing mesh data
	 */
	virtual void load_mesh(const std::string& filename) =0;

	///////////// Mesh properties

	/**
     * This method must always be implemented. The points may or may not
     * be duplicated, but all points in the mesh must be in the returned
     * container.
	 * @return set of active points/vertices of the mesh
	 */
	virtual const PointContainer* points() const =0;

	/**
     * This method is undefined if has_faces() returns false.
	 * @return set of vertex normals for each active point/vertex
	 * of the mesh. Let points be the vector returned from calling
	 * points() and normals be the vector returned from calling
	 * this method. Then normals[i] must be the vertex normal of
	 * the vertex at points[i]
	 */
	virtual const VertexNormalContainer* vertex_normals() const =0;

	/**
	 * @return set of active faces of the mesh. Let result be the
	 * returned vector. Then result[3*i], result[3*i+1] and
	 * result[3*i+2] contain the indices of the 3 vertices of face
	 * i in the vector returned by points().
	 */
	virtual const FaceContainer* faces() const =0;

    /**
     * @return true iff the implementation uses vertex indices.
     * i.e. is able to return the array of faces indexing the vertex
     * array returned by points(). Calls to points(), vertex_normals()
     * and faces() are valid if this is true
     */
    virtual bool has_vertex_indices() const =0;

    /**
     * Similar to the container returned by faces() but instead of
     * indexing the vertex array returned by points(), each element
     * references the vertex itself. This is for compatibility with
     * glDrawArrays call.
     * If the underlying implementation does not support this call,
     * it should return false in has_face_vertices(), in which case,
     * the pointer returned by this function is undefined.
     * @return set of vertices making up the faces of the mesh.
     */
    virtual const PointContainer* face_vertices() const =0;

    /**
     * This method is undefined if has_face_vertices() returns false.
     * @return set of face normals for each vertex returned by
     * face_vertices()
     */
    virtual const VertexNormalContainer* face_normals() const =0;

    /**
     * @return true iff the implementation supports the face_vertices()
     * and face_normals() call.
     */
    virtual bool has_face_vertices() const =0;

    /**
     * @return number of faces the mesh has
     */
    virtual size_t num_faces() const =0;
};

}

#endif
