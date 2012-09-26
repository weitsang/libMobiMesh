/*
 *  ProgressiveMesh.h
 *  MobiMesh
 *
 *  Created by Ong Yuh Shin on 3/9/11.
 *  Copyright 2011 National University of Singapore. All rights reserved.
 *
 */

#ifndef __PROGRESSIVE_MESH_H_
#define __PROGRESSIVE_MESH_H_

#include <string>
#include <stdint.h>

#include "IMesh.h"

namespace MobiMesh {

/**
 * Represents a progressive mesh object.
 * @see Progressive Meshes, Hoppe, SIGGRAPH 96
 */
class ProgressiveMesh : public IMesh
{
public:

	///////////// Constructors
	ProgressiveMesh() : n_base_vertices_(0), n_base_faces_(0), n_details_(0) { }
	virtual ~ProgressiveMesh() { }

	/**
	 * Loads the supplied file into the the object.
	 * @param[in] filename name of file containing pm data
	 * @throw InvalidFileFormatException if file referred to by
     * filename does not contain a valid PM data file
	 */
	explicit ProgressiveMesh(const std::string& filename);

	///////////// Methods inherited from interface IMesh

    /**
     * @see IMesh::load_mesh(const std::string& filename)
	 * @throw InvalidFileFormatException if file referred to by
     * filename does not contain a valid PM data file
     */
	virtual void load_mesh(const std::string& filename);

	const IMesh::PointContainer* points() const { return &points_; }

	const IMesh::VertexNormalContainer* vertex_normals() const { return &vertex_normal_; }

	const IMesh::FaceContainer* faces() const { return &draw_faces_; }

    bool has_vertex_indices() const { return true; }

    IMesh::PointContainer* face_vertices() const { return NULL; }

    IMesh::VertexNormalContainer* face_normals() const { return NULL; }

    bool has_face_vertices() const { return false; }

    size_t num_faces() const { return faces_.size(); }

	///////////// ProgressiveMesh specific methods

    virtual void unload_mesh();

protected:

	///////////// Halfedge data structures

	typedef int VertexHandle;
	typedef int HalfedgeHandle;
	typedef int EdgeHandle;
	typedef int FaceHandle;

	typedef struct Vertex
	{
		HalfedgeHandle halfedge_handle_;
		Vertex() : halfedge_handle_(-1) { }
	} Vertex;

	typedef struct Halfedge
	{
		FaceHandle     face_handle_;
		VertexHandle   vertex_handle_;
		HalfedgeHandle next_halfedge_handle_;
		HalfedgeHandle prev_halfedge_handle_;
		Halfedge() :
		  face_handle_(-1),
		  vertex_handle_(-1),
		  next_halfedge_handle_(-1),
		  prev_halfedge_handle_(-1) { }
	} Halfedge;

	typedef struct Edge
	{
		Halfedge halfedges_[2];
	} Edge;

	typedef struct Face
	{
		HalfedgeHandle halfedge_handle_;
		Face() : halfedge_handle_(-1) { }
	} Face;

	typedef std::vector<Vertex> VertexList;
	typedef std::vector<Edge>   EdgeList;
	typedef std::vector<Face>   FaceList;

	///////////// Adding vertices, faces, edges to the mesh

	/**
	 * Creates a new vertex at point (0,0,0).
	 *
	 * @return handle to the newly created vertex.
	 * @note subclasses may override this to provide memory management for
	 * vertices in the mesh.
	 */
	virtual VertexHandle new_vertex();

	/**
	 * Creates a new vertex at the given point.
	 *
	 * @param[in] p position of the new vertex
	 * @return handle to the newly created vertex
	 */
	virtual VertexHandle new_vertex(const IMesh::Point& p);

	/**
	 * Creates a new edge (with default start/end vertices in Edge's
	 * default constructor)
	 *
	 * @return halfedge handle of the halfedge in the direction from the
	 * start to the end vertex.
	 * @note subclasses may override this to provide memory management for
	 * edges in the mesh.
	 */
	virtual HalfedgeHandle new_edge();

	/**
	 * Creates a new edge connecting the two vertices supplied.
	 *
	 * @param[in] start_vertex_handle a valid vertex handle of the start
	 * vertex of the edge.
	 * @param[in] end_vertex_handle a valid vertex handle of the end vertex
	 * of the edge.
	 * @return halfedge handle of the halfedge going from start to end.
	 */
	virtual HalfedgeHandle new_edge(
	  VertexHandle start_vertex_handle,
	  VertexHandle end_vertex_handle);

	template <class _Handle>
		struct NextCacheEntryT : public std::pair<_Handle, _Handle>
	{
		typedef std::pair<_Handle, _Handle> Base;
		NextCacheEntryT(_Handle _heh0, _Handle _heh1) : Base(_heh0, _heh1)
		{ assert(_heh0 != -1); assert(_heh1 != -1); }
	};

	/**
	 * Creates a new face (with params using Face's default construcutor)
	 *
	 * @return face handle of the newly created face.
	 * @note subclasses may override this to provide memory management for
	 * edges in the mesh.
	 */
	virtual FaceHandle new_face();

	/**
	 * Creates a new face.
	 *
	 * @param[in] i0 a valid vertex handle to first vertex of the face.
	 * @param[in] i1 a valid vertex handle to second vertex of the face.
	 * @param[in] i3 a valid vertex handle to third vertex of the face.
	 * @return handle to newly created face, invalid handle if unable to
	 * create the face.
	 */
	virtual FaceHandle new_face(
	  VertexHandle i0, VertexHandle i1, VertexHandle i2);

	///////////// Accessors - Retrieves a reference to the object referenced
	///////////// by the supplied handle. Handle must be valid (otherwise
	///////////// we will access invalid memory).

	const Vertex& get_vertex(VertexHandle vertex_handle) const
	{ return      vertices_[vertex_handle]; }
	const OpenMesh::Vec3f& get_vertex_normal(VertexHandle vertex_handle) const
	{ return      vertex_normal_[vertex_handle]; }
	const Edge&   get_edge(EdgeHandle edge_handle) const
	{ return      edges_[edge_handle]; }
	const Halfedge& get_halfedge(HalfedgeHandle halfedge_handle) const
	{ return      edges_[halfedge_handle>>1].halfedges_[halfedge_handle&1]; }
	const Halfedge& get_opposite_halfedge(HalfedgeHandle halfedge_handle) const
	{ return      edges_[halfedge_handle>>1].halfedges_[(halfedge_handle&1)^1]; }
	const Face&   get_face(FaceHandle face_handle) const
	{ return      faces_[face_handle]; }
	const IMesh::Point& get_point(VertexHandle vertex_handle) const
	{ return      points_[vertex_handle]; }

	Vertex&   get_vertex(VertexHandle vertex_handle)
	{ return  vertices_[vertex_handle]; }
	OpenMesh::Vec3f& get_vertex_normal(VertexHandle vertex_handle)
	{ return  vertex_normal_[vertex_handle]; }
	Edge&     get_edge(EdgeHandle edge_handle)
	{ return  edges_[edge_handle]; }
	Halfedge& get_halfedge(HalfedgeHandle halfedge_handle)
	{ return  edges_[halfedge_handle>>1].halfedges_[halfedge_handle&1]; }
	Halfedge& get_opposite_halfedge(HalfedgeHandle halfedge_handle)
	{ return      edges_[halfedge_handle>>1].halfedges_[(halfedge_handle&1)^1]; }
	Face&     get_face(FaceHandle face_handle)
	{ return  faces_[face_handle]; }
	IMesh::Point& get_point(VertexHandle vertex_handle)
	{ return  points_[vertex_handle]; }

	///////////// Checks to see if a handle refers to a valid object.

	bool     is_valid_vertex(VertexHandle vertex_handle) const
	{ return vertex_handle != -1; }
	bool     is_valid_edge(EdgeHandle edge_handle) const
	{ return edge_handle != -1; }
	bool     is_valid_halfedge(HalfedgeHandle halfedge_handle) const
	{ return halfedge_handle != -1; }
	bool     is_valid_face(FaceHandle face_handle) const
	{ return face_handle != -1; }

	///////////// Geometry-related accessors - convenient methods to
	///////////// access and manipulate vertices/edges/faces/etc.

	/**
	 * @param[in] heh halfedge handle of the halfedge whose face handle
	 * we want. Must be a valid handle.
	 * @return face handle of the face the halfedge heh references.
	 */
	FaceHandle face_handle(HalfedgeHandle heh) const
	{ return get_halfedge(heh).face_handle_; }

	/**
	 * @param[in] heh halfedge handle of the edge whose handle we want.
	 * Must be a valid handle.
	 * @return edge handle of the edge corresponding to heh
	 */
	EdgeHandle edge_handle(HalfedgeHandle heh) const
	{ return heh >> 1; }

	/**
	 * @param[in] heh halfedge handle of the halfedge whose face handle
	 * we want to change. Must be a valid handle.
	 * @param[in] fh face handle of the face the halfedge heh should
	 * reference.
	 */
	void set_face_handle(HalfedgeHandle heh, FaceHandle fh)
	{ get_halfedge(heh).face_handle_ = fh; }

	/**
	 * @param[in] vh vertex handle of the vertex whose halfedge we want.
	 * Must be a valid handle.
	 * @return halfedge handle of the halfedge referenced by this vh.
	 */
	HalfedgeHandle halfedge_handle_vh(VertexHandle vh) const
	{ return get_vertex(vh).halfedge_handle_; }

	/**
	 * @param[in] vh vertex handle of the vertex whose halfedge we want
	 * to change. Must be a valid handle.
	 * @param[in] heh halfedge handle of the halfedge that vertex vh should
	 * reference.
	 */
	void set_halfedge_handle_vh(VertexHandle vh, HalfedgeHandle heh)
	{ get_vertex(vh).halfedge_handle_ = heh; }

	/**
	 * @param[in] fh face handle of the face whose halfedge we want.
	 * Must be a valid handle.
	 * @return halfedge handle of the halfedge referenced by fh.
	 */
	HalfedgeHandle halfedge_handle_fh(FaceHandle fh) const
	{ return get_face(fh).halfedge_handle_; }

	/**
	 * @param[in] fh face handle of the face whose halfedge we want
	 * to change. Must be a valid handle.
	 * @param[in] heh halfedge handle of the halfedge that fh should
	 * reference.
	 */
	void set_halfedge_handle_fh(FaceHandle fh, HalfedgeHandle heh)
	{ get_face(fh).halfedge_handle_ = heh; }

	/**
	 * @param[in] heh halfedge handle of the halfedge whose next halfedge
	 * we want. Must be a valid handle.
	 * @return halfedge handle to the next halfedge pointed to by heh.
	 */
	HalfedgeHandle next_halfedge_handle(HalfedgeHandle heh) const
	{ return get_halfedge(heh).next_halfedge_handle_; }

	/**
	 * @param[in] heh halfedge handle of the halfedge whose previous
	 * halfedge we want. Must be a valid handle.
	 * @return halfedge handle to the previous halfedge pointed to by heh.
	 */
	HalfedgeHandle prev_halfedge_handle(HalfedgeHandle heh) const
	{ return get_halfedge(heh).prev_halfedge_handle_; }

	/**
	 * Links up heh to nheh, modifying their next and previous halfedge
	 * references respectively.
	 * @param[in] heh halfedge handle of the halfedge whose next halfedge
	 * is to be changed. Must be a valid handle.
	 * @param[in] nheh halfedge handle of next halfedge.
	 */
	void set_next_halfedge_handle(
		HalfedgeHandle heh, HalfedgeHandle nheh)
	{
		get_halfedge(heh).next_halfedge_handle_ = nheh;
		get_halfedge(nheh).prev_halfedge_handle_ =  heh;
	}

	/**
	 * @param[in] halfedge_handle halfedge handle to the halfedge whose
	 * opposite we want to find. Must be a valid handle.
	 * @return the halfedge handle of the halfedge opposite to the
	 * given one
	 */
	HalfedgeHandle get_opposite_halfedge_handle(
	  HalfedgeHandle halfedge_handle) const
	{ return (halfedge_handle & 1) ? halfedge_handle-1 : halfedge_handle+1; }

	/**
	 * @param[in] start_vertex_handle vertex handle of the vertex where
	 * the halfedge originates from. Must be a valid handle.
	 * @param[in] end_vertex_handle vertex handle of the vertex where the
	 * halfedge terminates. Must be a valid handle.
	 * @return halfedge handle for halfedge from start vertex to end vertex,
	 * or an invalid handle if no such edge exists
	 */
	HalfedgeHandle find_halfedge(
	  VertexHandle start_vertex_handle,
	  VertexHandle end_vertex_handle) const;

	/**
	 * A boundary halfedge is one which is not adjacent to any face.
	 * Note that it's opposite halfedge can be adjacent to another face.
	 * @param[in] halfedge_handle handle to the halfedge to be checked,
	 * must be a valid handle.
	 */
	bool is_boundary_halfedge(HalfedgeHandle halfedge_handle) const
	{ return !is_valid_face(get_halfedge(halfedge_handle).face_handle_); }

	/**
	 * A boundary vertex is one which is not associated with a valid
	 * halfedge or face.
	 * @param[in] vertex_handle handle to the vertex to be checked,
	 * must be a valid handle.
	 */
	bool is_boundary_vertex(VertexHandle vertex_handle) const
	{
		HalfedgeHandle heh = vertices_[vertex_handle].halfedge_handle_;
		if (!is_valid_halfedge(heh)) return true;
		return !is_valid_face(get_halfedge(heh).face_handle_);
	}

	///////////// Get containers' sizes
    
	size_t n_faces()     const { return faces_.size(); }
	size_t n_vertices()  const { return vertices_.size(); }
	size_t n_halfedges() const { return  2*edges_.size(); }
	size_t n_edges()     const { return    edges_.size(); }

	///////////// PM vertex split operation

	/**
	 * Performs a vertex split operation, typically introducing 2 new
	 * faces, unless splitting a boundary vertex, where either vl or vr
	 * is invalid, in which case, only 1 face is introduced.
	 * Halfedges are added and re-assigned accordingly.
	 *
	 * @note vl and vr cannot both be invalid
	 * @note the caller is responsible for creating the new vertex v0
	 * (and passing in a valid handle)
	 *
	 * @param[in] v0 vertex handle of newly added vertex due to split.
	 * @param[in] v1 vertex handle of vertex to be be split,
	 * must not be invalid
	 * @param[in] vl vertex handle of vertex that will be opposite v0v1
	 * in the newly added left face, invalid handle if vl does not exist
	 * @param[in] vr vertex handle of vertex that will be opposite v0v1
	 * in the newly added right face, invalid handle if vr does not exist
	 * @param[in] v1vl halfedge handle of the halfedge from v1 to vl
	 * @param[in] vrv1 halfedge handle of the halfedge from vr to v1
	 * @return halfedge handle to newly added edge in direction v0->v1
	 */
	HalfedgeHandle vertex_split(
	  VertexHandle v0, VertexHandle v1, VertexHandle vl, VertexHandle vr,
		HalfedgeHandle v1vl, HalfedgeHandle vrv1);

    ///////////// The following should be private data members. However,
    ///////////// for ease of implementation in child classes (e.g.
    ///////////// garbage collection and edge collapse operations), these
    ///////////// are exposed. Where possible, accessors to specific
    ///////////// geometry data should be made via the get_* methods.

	/// GPU friendly face array
	IMesh::FaceContainer  draw_faces_;
    
	IMesh::PointContainer        points_;
	IMesh::VertexNormalContainer vertex_normal_;

    VertexList vertices_;
	EdgeList   edges_;
	FaceList   faces_;

private:

	///////////// Copying disallowed (enable only if really needed, use
	///////////// pointers or references otherwise)
	ProgressiveMesh(const ProgressiveMesh&);
	ProgressiveMesh& operator=(const ProgressiveMesh&);

	///////////// PM implementation specific methods and data members

	/**
	 * Recompute normals at each vertex of each face
	 */
	void update_normals();

	/**
	 * Finds a boundary halfedge that originates from the vertex referred
	 * to by vertex_handle (if any) and sets the vertex's halfedge handle
	 * to be the boundary halfedge. Nothing is changed if the vertex is
	 * not adjacent to any such boundary half edge. This is so that
	 * boundary vertices can be detected by simply checking whether their
	 * halfedge is a boundary edge.
	 *
	 * @param[in] vertex_handle vertex handle of the vertex whose halfedge
	 * handle needs to be set to some boundary halfedge (if any)
	 */
	void adjust_outgoing_halfedge(VertexHandle vertex_handle);

	/**
	 * Builds a loop between 2 connected vertices vl and v1: vl->v1->vl.
	 * @param[in] v1 vertex handle of v1.
	 * @param[in] vl vertex handle of vl.
	 * @param[in] v1vl halfedge handle of the existing edge from v1 to vl.
	 * @return halfedge handle of the newly created edge in the loop, from
	 * vl to v1.
	 */
	HalfedgeHandle build_loop(
		VertexHandle v1, VertexHandle vl, HalfedgeHandle v1vl);

	/// Number of vertices in base mesh
	uint32_t n_base_vertices_;
	/// Number of faces in base mesh
	uint32_t n_base_faces_;
	/// Number of details (vertex splits)
	uint32_t n_details_;
};

}

#endif
