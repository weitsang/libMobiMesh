/*
 *  ViewDependentProgressiveMesh.h
 *  MobiMesh
 *
 *  Created by Ong Yuh Shin on 3/10/11.
 *  Copyright 2011 National University of Singapore. All rights reserved.
 *
 */

#ifndef __VIEW_DEPENDENT_PROGRESSIVE_MESH_H_
#define __VIEW_DEPENDENT_PROGRESSIVE_MESH_H_

#include <list>

#include <unistd.h>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/mman.h>
#include <sys/stat.h>

#include "ProgressiveMesh.h"
#include "StatusInfo.h"
#include "ViewingParameters.h"

namespace MobiMesh {

/**
 * Represents a view-dependent progressive mesh (VPM).
 * @see View-dependent Refinement of Progressive Meshes,
 * Hoppe, SIGGRAPH 97
 * @see Truly Selective Refinement of Progressive Meshes,
 * Kim, Lee, GRIN'01
 */
class ViewDependentProgressiveMesh : public ProgressiveMesh
{
public:

	///////////// Constructors
	ViewDependentProgressiveMesh();
	virtual ~ViewDependentProgressiveMesh();

	/**
	 * Loads the supplied file into the object
	 * @param[in] filename name of file containing vpm data
	 * @throw FileReadErrorException if unable to read file
	 * @throw InvalidFileFormatException if file does not refer to a
     * VPM data file
	 */
	explicit ViewDependentProgressiveMesh(const std::string& filename);

	///////////// Methods overriding those in ProgressiveMesh

	virtual void load_mesh(const std::string& filename);
    void unload_mesh();
    bool is_mesh_loaded() { return is_mesh_loaded_; }

	///////////// Methods specific to VPM

	/**
	 * Perform adaptive refinement on the mesh based on current viewing
     * parameters.
     * @param max_refinement max number of refinements allowed, i.e. no more
     * than max_refinement vertex splits, and max_refinement edge collapses
     * will take place.
	 * @return true iff the mesh (i.e. faces, vertices, etc.) was refined.
	 */
	bool adapt(int max_refinement);

	/**
	 * Can be used to perform a check before calling adapt() so we do not call adapt
	 * unnecessarily.
	 * @return true iff the mesh needs to be refined.
	 */
	bool need_refined() { return need_refined_; }

	void set_viewing_parameters(const ViewingParameters& viewing_parameters)
	{ viewing_parameters_ = viewing_parameters; }

    /**
     * Puts the mesh in active state and loads the necessary resources
     * so that adapt() can be called. If mesh is inactive, adapt() cannot
     * be called. Geometry accessor methods are still valid though because
     * the vertex front is preserved even when deactivated.
     */
    void activate();

    /**
     * Unloads any resources that will not be needed by the mesh if
     * it does not need to be refined. In this state, the behavior of
     * adapt() is undefined. However, the current vertex front is still
     * active and valid if the caller attempts to access them.
     *
     * The difference between activating and loading a mesh is that
     * activation/deactivation preserves the current vertex front while
     * unloading clears the internal structures entirely. Reloading the
     * mesh after unloading takes a long time, whereas re-activating
     * should be quick.
     */
    void deactivate();

    bool is_active() { return is_active_; }

protected:

	///////////// Vertex parameters (for nodes in vpm's vertex tree)

	typedef unsigned int VertexID;

	/**
	 * Vertex parameters for each vertex in the vertex tree.
	 */
	typedef struct VParam
	{
		/// Identifier for node
		VertexID vid;

		/// Handle to the actual vertex object
		VertexHandle vh;

		/// Index of vertex split info for this vertex
		unsigned pos;

		/// Vertex normal
		OpenMesh::Vec3f normal;

		/// View-dependent params
		float radius, sin_square, mue_square, sigma_square;

		/// Pointer to parent vertex in vertex tree
		std::list<VParam>::iterator parent;

		/// Pointer to active children (if any) in vertex tree
		std::list<VParam>::iterator lchild, rchild;
	} VParam;

	/**
	 * Connectivity & View-dependent parameters for a vertex split.
	 * This struct represents how the VPM is stored in the file,
	 * i.e. if offset o points to the start of a vertex split data in
	 * the file, then the next o + sizeof(VsData) bytes can be read
	 * directly into an instance of this struct to extract the vertex
	 * split information.
	 */
	typedef struct VsData {
		/// Position of vertex to be added due to this split
		OpenMesh::Vec3f p;

		/// Vertex ID (in vertex tree) of this vertex
		unsigned int node_index;

		/// Vertex ID (in vertex tree) of left cut neighbour
		unsigned int fund_lcut_index;

		/// Vertex ID (in vertex tree) of right cut neighbour
		unsigned int fund_rcut_index;

		// Vertex params for the left vertex from the split
		float l_radius;
		OpenMesh::Vec3f l_normal;
		float l_sin_square;
		float l_mue_square;
		float l_sigma_square;

		// Vertex params for the right vertex from the split
		float r_radius;
		OpenMesh::Vec3f r_normal;
		float r_sin_square;
		float r_mue_square;
		float r_sigma_square;

		/// Index of VsData for splitting the resulting left vertex
		unsigned int l_pos;

		/// Index of VsData for splitting the resulting right vertex
		unsigned int r_pos;
	} VsData;

	typedef int Remark;

	typedef std::vector<std::list<VParam>::iterator> VParamContainer;
	typedef std::vector<StatusInfo>                  VStatusContainer;
	typedef std::vector<StatusInfo>                  EStatusContainer;
	typedef std::vector<Normal>                      FNormalContainer;
	typedef std::vector<StatusInfo>                  FStatusContainer;
	typedef std::vector<Remark>                      VRemarkContainer;

	///////////// Overriding methods in PM class

	virtual VertexHandle new_vertex();
	virtual VertexHandle new_vertex(const Point& p);

	virtual HalfedgeHandle new_edge();
	virtual HalfedgeHandle new_edge(
	  VertexHandle start_vertex_handle,
	  VertexHandle end_vertex_handle);

	virtual FaceHandle new_face();
	virtual FaceHandle new_face(
	  VertexHandle i0, VertexHandle i1, VertexHandle i2);

	///////////// VPM specific methods

	/**
	 * Refine the mesh by recomputing the active faces and vertices of based on the
	 * current view. Caller needs to update the viewing parameters prior to calling
	 * this method to adapt the mesh to the current view.
	 * The number of edge collapses/vertex splits is limited in this call (so that
	 * the caller can possibly have a better frame rate)
	 *
	 * @return true iff the mesh will need further refinement
	 */
	bool adaptive_refinement(int max_refinement);

	/**
	 * Performs a vertex split.
	 *
	 * @param[in] vpit iterator pointing to the vertex param object
	 * describing the vertex to be split.
     * @return true iff split was performed. Vertex split will not be done if
     * there are too many faces or vertices
	 */
	bool force_vsplit(std::list<VParam>::iterator vpit);
	
	/**
	 * Performs an edge collapse operation. Updates and removes the
	 * vertex and its params accordingly
	 *
	 * @param[in] vpit the vertex param of the vertex whose vertex
	 * split introduces the edge that we want to collapse
	 * @param[in] v0v1 halfedge handle of the edge to be collapsed
	 */
	void ecol(std::list<VParam>::iterator vpit, HalfedgeHandle v0v1);

	///////////// Container accessor methods - returns an iterator/
	///////////// pointer to the object in each container based on the
	///////////// supplied handle. Standard * and -> operators apply
	///////////// if the pointer/iterator is valid. Pointer/iterator
	///////////// is valid iff supplied handle is valid.

	typedef VParamContainer::iterator  VParamPointer;
	typedef VStatusContainer::iterator VStatusPointer;
	typedef EStatusContainer::iterator EStatusPointer;
	typedef FNormalContainer::iterator FNormalPointer;
	typedef FStatusContainer::iterator FStatusPointer;

	VParamPointer  get_vertex_param(VertexHandle vertex_handle) const;
	VStatusPointer get_vertex_status(VertexHandle vertex_handle) const;
	EStatusPointer get_edge_status(EdgeHandle edge_handle) const;
	FNormalPointer get_face_normal(FaceHandle face_handle) const;
	FStatusPointer get_face_status(FaceHandle face_handle) const;

private:

	///////////// Copying disallowed (enable only if really needed, use
	///////////// pointers or references otherwise)
	ViewDependentProgressiveMesh(const ViewDependentProgressiveMesh&);
	ViewDependentProgressiveMesh& operator=(const ViewDependentProgressiveMesh&);
	
	///////////// VPM implementation specific methods and data members

	/**
	 * Compute and store the number of tree id bits required to represent the given
	 * number of trees. Assuming a node id is 32-bit.
	 */
	void calc_tree_id_bits(unsigned int num_roots);

	/**
	 * Queries whether a vertex (with params vp) should be split based
	 * on the current view parameters.
	 *
	 * @param[in] vp vertex param of the vertex to be split
	 * @return true iff the vertex identified by vp should be split
	 */
	bool qrefine(const VParam& vp);

	/**
	 * Checks if a vertex v's normal is oriented away from the viewer,
	 * if it is, it means the face associated with this vertex is not
	 * facing the viewer directly.
	 *
	 * @param[in] sin_square square of alpha_v, where alpha_v is the
	 * semiangle of a cone about v's normal
	 * @param[in] distance_square || v - e || ^ 2, where v is the
	 * position vector of the vertex, and e is the position vector of
	 * the eye
	 * @param[in] product_value (v - e).n_v, where n_v is the unit
	 * normal vector of vertex v
	 */
	bool oriented_away(
	  float sin_square, float distance_square, float product_value) const;

	/**
	 * Measures whether the distance between the approximate surface
	 * and the surface in the original mesh, when projected on the screen,
	 * is everywhere less than a screen space error tolerance.
	 *
	 * @param[in] mue_square square of the uniform component of the
	 * deviation space
	 * @param[in] sigma_square square of the radius of the deviation space
	 * @param[in] distance_square || v - e || ^2 where v is the vertex
	 * position and e is the eye position
	 * @param[in] product_value (v-e).n_v where n_v is the unit normal
	 * vector for the vertex at v
	 * @return true iff the screen space projection error exceeds the
	 * error tolerance constant
	 */
	bool screen_space_error(
	  float mue_square, float sigma_square, 
	  float distance_square, float product_value) const;

	/**
	 * Returns left and right cut of vertex referred to by vpit.
	 * The fundamental/actual left and right cut vertices might not
	 * be active. In which case, we return the active vertex that is
	 * an ancestor of the left/right cut vertex respectively.
	 *
	 * @param[in] vpit iterator to the vertex parameter of the vertex
	 * whose left and right cut handles we want to find
	 * @param[in] vs vertex split data for the vertex *vpit
	 * @param[out] vl the vertex handle to the vertex that is either
	 * the left cut of the vertex referred to by vpit (if the left cut
	 * exists in the current mesh), or the vertex that needs to be split
	 * eventually in order to obtain this left cut
	 * @param[out] vr similar to vl but for the right cut
	 * @param[out] v1vl halfedge handle of halfedge from v1 to vl
	 * @param[out] vrv1 halfedge handle of halfedge from vr to v1
	 */
	void get_active_cuts(std::list<VParam>::iterator vpit, const VsData& vs,
                         VertexHandle& vl, VertexHandle& vr,
                         HalfedgeHandle& v1vl, HalfedgeHandle& vrv1) const;

	/**
	 * Performs a vertex split. The necessary faces and halfedges will
	 * be added/re-assigned accordingly. Vertex params are created and
	 * updated accordingly.
	 *
	 * @param[in,out] vpit iterator to the object containing the vertex
	 * params of the vertex to be split. Upon returning vpit will have
	 * its lchild and rchild values updated to contain the VParam info of
	 * the newly added vertices.
	 * @param[in] vs vertex split info
	 * @param[in] vertex handle to the left cut vertex of the vertex split
	 * @param[in] vertex handle to the right cut vertex of the vertex split
     * @return true iff split was performed. Vertex split will not be done if
     * there are too many faces or vertices
	 */
	bool vsplit(std::list<VParam>::iterator vpit, const VsData& vs,
                VertexHandle vl, VertexHandle vr,
                HalfedgeHandle v1vl, HalfedgeHandle vrv1);

	/**
	 * Checks if it is ok to collapse the edge v0v1.
	 *
	 * @param[in] v0 handle of the vertex v0
	 * @param[in] v1 handle of the vertex v1
	 * @param[in] v0v1 handle to half-edge from vertex v0 to vertex v1,
	 * where v0v1 is the edge to be collapsed
	 * @return true iff the collapsing of this edge is possible. Sometimes
	 * an edge cannot be collapsed because collapsing it causes mesh problems
	 * such as a face being made up by only 2 vertices.
	 */
	bool is_collapse_ok(VertexHandle v0, VertexHandle v1, HalfedgeHandle& v0v1);

	/**
	 * @param hh halfedge handle of the edge to be collapsed
	 */
	void collapse(HalfedgeHandle hh);

	/**
	 * Removes a loop comprised of halfedges hp -> hn -> hp
	 */
	void collapse_loop(
		HalfedgeHandle hp, Halfedge& hpe, HalfedgeHandle hn, Halfedge& hne);

	/**
	 * Removes vertices, edges, faces that are flagged for deletion.
	 * @param[in] v true iff vertices flagged for deletion are to be cleaned up.
	 * @param[in] e true iff edges flagged for deletion are to be cleaned up.
	 * @param[in] f true iff faces flagged for deletion are to be cleaned up.
	 */
	void garbage_collection(bool v, bool e, bool f);

	/**
	 * Loads vertex split data at the given index.
	 */
	VsData loadVs(unsigned int data_index);

	bool is_ancestor_VertexID(VertexID ancestor, VertexID descendent) const;

	VParamContainer   vertex_params_;
	VStatusContainer  vertex_status_;
	VRemarkContainer  vertex_remark_;
	EStatusContainer  edge_status_;
	FNormalContainer  face_normal_;
	FStatusContainer  face_status_;

	ViewingParameters viewing_parameters_;
	float             kappa_square_;

	// number of bits used to represent the tree id
	unsigned char tree_id_bits_;
	unsigned char tree_id_bitshift_; //          32 - tree_id_bits_
	unsigned int  tree_id_mask_;     // 0xFFFFFFFF >> tree_id_bits_

	// use mmap to map files and manage caches
	int         fd_;
	char*       data_;
	struct stat sbuf_;

    /// whether mmap() call has been made and has succeeded
    bool mmapped_;

	/// number of bytes the base mesh information takes (in the .vpm file)
	size_t base_offset_;
	/// number of bytes each vertex split data takes (in the .vpm file)
	size_t data_size_;
 
	typedef std::list<VParam>::iterator VPIter;

	/// maximum number of edge collapse/vertex splits per call to adaptive refinement
	static const unsigned int MAX_ECVS_PER_REFINEMENT = 2000;

	/// list of free/reusable vertices
	std::vector<size_t> free_vertices_;
	/// smallest index of a reusable vertex in free_vert_
	int free_vertex_index_;

	/// list of free/reusable faces
	std::vector<size_t> free_faces_;
	/// smallest index of a reusable face in free_faces_
	int free_face_index_;

	/// list of free/reusable edges
	std::vector<size_t> free_edges_;
	/// smallest index of a reusable edge in free_edges_
	int free_edge_index_;

	/// list of free/reusable vertex params
	std::vector<VPIter> free_vplist_;
	/// smallest index of a reusable vertex param in free_vplist_
	int free_param_index_;

	/// list of vertex params for all vertices in the active vertex tree
	std::list<VParam> vplist_;

	/// flag to determine whether a mesh is currently loaded
	bool is_mesh_loaded_;

    /// flag determining whether mesh is active, i.e. if adapt() can be called
    bool is_active_;

    /// filename of the file containing the mesh
    std::string mesh_filename_;

	/// Number of vertices in base mesh
	uint32_t n_base_vertices_;
	/// Number of faces in base mesh
	uint32_t n_base_faces_;
	/// Number of details (vertex splits)
	uint32_t n_details_;

	/// index/vertex handle of the last vertex that was last processed for split/collapse
	size_t vit_;

	/// whether the mesh needs to be further refined based on the last refinement and
	/// view parameters.
	bool need_refined_;
};

}

#endif
