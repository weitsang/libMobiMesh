/*
 *  ViewDependentProgressiveMesh.cpp
 *  MobiMesh
 *
 *  Created by Ong Yuh Shin on 3/10/11.
 *  Copyright 2011 National University of Singapore. All rights reserved.
 *
 */

#include <fstream>
#include <fcntl.h>
#include <errno.h>
#include <unistd.h>

#include "MobiMeshConstants.h"
#include "FileReadErrorException.h"
#include "InvalidFileFormatException.h"
#include "VHierarchyNodeIndex.h"
#include "ViewDependentProgressiveMesh.h"

namespace MobiMesh {

ViewDependentProgressiveMesh::ViewDependentProgressiveMesh() :
	kappa_square_(0.0), tree_id_bits_(0), fd_(0), data_(NULL),
	base_offset_(0), data_size_(0),
	free_vertex_index_(-1), free_face_index_(-1),
	free_edge_index_(-1), free_param_index_(-1),
	is_mesh_loaded_(false), need_refined_(false),
    mmapped_(false), is_active_(false)
{
	free_vertices_.reserve(MAX_ECVS_PER_REFINEMENT);
	free_faces_.reserve(2*MAX_ECVS_PER_REFINEMENT);
	free_edges_.reserve(3*MAX_ECVS_PER_REFINEMENT);
	free_vplist_.reserve(2*MAX_ECVS_PER_REFINEMENT);
}

ViewDependentProgressiveMesh::ViewDependentProgressiveMesh(const std::string& filename)
	: kappa_square_(0.0), tree_id_bits_(0), fd_(0), data_(NULL),
	base_offset_(0), data_size_(0),
	free_vertex_index_(-1), free_face_index_(-1),
	free_edge_index_(-1), free_param_index_(-1),
	is_mesh_loaded_(false), need_refined_(false),
    mmapped_(false), is_active_(false)
{
	free_vertices_.reserve(MAX_ECVS_PER_REFINEMENT);
	free_faces_.reserve(2*MAX_ECVS_PER_REFINEMENT);
	free_edges_.reserve(3*MAX_ECVS_PER_REFINEMENT);
	free_vplist_.reserve(2*MAX_ECVS_PER_REFINEMENT);

	load_mesh(filename);
}

ViewDependentProgressiveMesh::~ViewDependentProgressiveMesh()
{
	if (is_mesh_loaded_)
		unload_mesh();
}

void ViewDependentProgressiveMesh::load_mesh(const std::string& filename)
{
	if (is_mesh_loaded_)
		unload_mesh();

	char            fileformat[16];
	unsigned int    fvi[3], pos;
	OpenMesh::Vec3f p, normal;
	float           radius, sin_square, mue_square, sigma_square;

	std::ifstream ifs(filename.c_str(), std::ios::binary);
	if (!ifs)
		throw FileReadErrorException("ifstream read error");

	ifs.read(fileformat, 10);
	fileformat[10] = '\0';
	if (std::string(fileformat) != std::string("VDProgMesh"))
	{
		ifs.close();
		throw InvalidFileFormatException(
			"Invalid file format", "VDProgMesh", fileformat);
	}

	ifs.read( (char*)&n_base_vertices_, sizeof(n_base_vertices_) );
	ifs.read( (char*)&n_base_faces_   , sizeof(n_base_faces_) );
	ifs.read( (char*)&n_details_      , sizeof(n_details_) );

	calc_tree_id_bits(n_base_vertices_);

	for (size_t i=0; i<n_base_vertices_; ++i)
	{
		ifs.read( (char*)&p            , sizeof(p) );
		ifs.read( (char*)&radius       , sizeof(radius) );
		ifs.read( (char*)&normal       , sizeof(normal) );
		ifs.read( (char*)&sin_square   , sizeof(sin_square) );
		ifs.read( (char*)&mue_square   , sizeof(mue_square) );
		ifs.read( (char*)&sigma_square , sizeof(sigma_square) );
		ifs.read( (char*)&pos          , sizeof(pos) );

		VParam vp;
		vp.vh            = new_vertex(p);
		vp.pos           = pos;
		vp.vid           = VHierarchyNodeIndex(i, 1, tree_id_bits_).value();
		vp.normal        = normal;
		vp.radius        = radius;
		vp.sin_square    = sin_square;
		vp.mue_square    = mue_square;
		vp.sigma_square  = sigma_square;
		vp.parent        = vplist_.end();
		vp.lchild        = vplist_.end();
		vp.rchild        = vplist_.end();
		
		get_vertex_normal(vp.vh) = vp.normal;
		vertex_params_[vp.vh] = vplist_.insert(vplist_.end(),vp);
	}

	for (size_t i=0; i<n_base_faces_; ++i)
	{
		ifs.read( (char*)&fvi[0] , sizeof(fvi[0]) );
		ifs.read( (char*)&fvi[1] , sizeof(fvi[1]) );
		ifs.read( (char*)&fvi[2] , sizeof(fvi[2]) );
		new_face(fvi[0],fvi[1],fvi[2]);
	}

	ifs.close();
	
	// We use vertex normals from .vpm directly without real-time updates; (0.5.2)
	// update_face_normals();
	// We calculate the base offset and data size to access mapped files properly; (0.5.7)
	size_t s_params_, s_base_vertices_, s_base_faces_;
	s_params_        = sizeof("VDProgMesh") - 1 + // size of mesh type in header
	                   sizeof(unsigned int) +  // size to store num vertices in base mesh
	                   sizeof(unsigned int) + // size to store num faces in base mesh
	                   sizeof(unsigned int); // size to store num vertex splits in vpm
	s_base_vertices_ = sizeof(OpenMesh::Vec3f) + // 3D coords of vertex
	                   sizeof(float) + // radius of bounding sphere
	                   sizeof(OpenMesh::Vec3f) + // 3D vertex normal
	                   sizeof(float) + // sin^2 alpha
	                   sizeof(float) + // mu^2
	                   sizeof(float) + // delta^2
	                   sizeof(unsigned int); // VsData index
	s_base_faces_    = sizeof(unsigned int)+sizeof(unsigned int)+sizeof(unsigned int);

	base_offset_     = s_params_+n_base_vertices_*s_base_vertices_ +
	                   n_base_faces_*s_base_faces_;
	data_size_       = sizeof(OpenMesh::Vec3f) +
	                   sizeof(unsigned int) +
	                   sizeof(unsigned int) +
	                   sizeof(unsigned int) +
	                   sizeof(float) +
	                   sizeof(OpenMesh::Vec3f) +
	                   sizeof(float) +
	                   sizeof(float) +
	                   sizeof(float) +
	                   sizeof(float) +
	                   sizeof(OpenMesh::Vec3f) +
	                   sizeof(float) +
	                   sizeof(float) +
	                   sizeof(float) +
	                   sizeof(unsigned int) +
	                   sizeof(unsigned int);

#ifdef DEBUG
	assert (data_size_ == sizeof(VsData));
#endif

    mesh_filename_ = filename;
    activate();

	is_mesh_loaded_ = true;
	need_refined_ = true;
}

void ViewDependentProgressiveMesh::activate()
{
    if (is_active())
        return;

    // Activating the mesh opens the file to be readable from disk
    // so that vertx splits can be loaded.

    if ((fd_ = open(mesh_filename_.c_str(), O_RDONLY)) == -1)
		throw FileReadErrorException("open() failed, unable to open file");
	if (stat(mesh_filename_.c_str(), &sbuf_) == -1)
    {
        close(fd_);
		throw FileReadErrorException("stat() read error");
    }
	if ((data_ = (char *) mmap(NULL, sbuf_.st_size, PROT_READ,
                               MAP_SHARED, fd_, 0)) != MAP_FAILED)
	{
#ifdef DEBUG
        std::cout << "mmap() succeeded" << std::endl;
#endif
        mmapped_ = true;
	}
    else
    {
#ifdef DEBUG
        std::cout << "mmap() failed, vertex splits will be read with "
                     "pread() instead" << std::endl;
#endif
        mmapped_ = false;
        data_ = NULL;
    }

    is_active_ = true;
}

void ViewDependentProgressiveMesh::deactivate()
{
    if (!is_active())
        return;

    is_active_ = false;

	if (mmapped_ && munmap(data_, sbuf_.st_size) == -1)
    {
#ifdef DEBUG
		std::cerr << "Error un-mmapping" << std::endl;
#endif
        mmapped_ = false;
    }

    data_ = NULL;

	close(fd_);
}

bool ViewDependentProgressiveMesh::adapt(int max_refinement)
{
    max_refinement = (max_refinement > MAX_ECVS_PER_REFINEMENT) ? MAX_ECVS_PER_REFINEMENT : max_refinement;
	need_refined_ = adaptive_refinement(max_refinement);

	return need_refined_;
}

bool ViewDependentProgressiveMesh::adaptive_refinement(int max_refinement)
{
	// switch the execution order of vertex-split and edge_collapse; (0.5.1)
	// use 'need_further_refined_' to reduce operations every frame and thus
	// improve framerate; (0.5.7)

	bool need_further_refined_ = false;

	// edge-collapse (not root)
	HalfedgeHandle v0v1 = -1;
	size_t eccount = 0;
	if (vit_ >= n_vertices())
		vit_ = 0;
	for ( ; vit_ != n_vertices(); )
	{
		// limit the number of edge-collapses every frame to improve the framerate; (0.5.7)
		if (eccount >= max_refinement)
		{
			need_further_refined_ = true;
			break;
		}

		std::list<VParam>::iterator  vpit = vertex_params_[vit_];
		std::list<VParam>::iterator pvpit = vpit->parent;

		if ( vpit != vplist_.end() && pvpit != vplist_.end())
		{
			pvpit->vh = pvpit->rchild->vh;
			// only query the parent from rchild (not lchild) to save some visual
			// queries; (0.6.1)
			if ((vpit == pvpit->rchild) // is_parent
				&& pvpit->lchild->lchild == vplist_.end()
				&& pvpit->rchild->lchild == vplist_.end()
				&& qrefine(*pvpit) != true
				&& is_collapse_ok(pvpit->lchild->vh, pvpit->rchild->vh, v0v1) == true)
			{
				ecol(pvpit, v0v1);
                ++eccount;
			}
			else
				++vit_;
		}
		else
			++vit_;
	}

	// update VertexHandle links since 'garbage_collection()' will shuffle them; (0.1.0)
	// move this routine into 'garbage_collection()' to be more efficient; (0.5.7)
	size_t vscount = 0;
	if (vit_ >= n_vertices())
		vit_ = 0;
	for ( ; vit_ != n_vertices(); )
	{
		if (n_vertices() >= Constants::MAX_VERTICES)
            break;

		// limit the number of vertex-splits every frame to improve the framerate; (0.5.7)
		if (vscount >= max_refinement)
		{
			need_further_refined_ = true;
			break;
		}

		std::list<VParam>::iterator vpit = vertex_params_[vit_];
		if ( vpit != vplist_.end() && vpit->pos != 0 && qrefine(*vpit))
        {
			if (force_vsplit(vpit))
                ++vscount;
        }
		else
			++vit_;
	}

	// move these resize routines here to speed up a single vertex-split; (0.1.0)
	edge_status_.resize(n_edges());
	face_normal_.resize(n_faces());
	face_status_.resize(n_faces());

	// collect those elements (v/e/f) tagged as 'deleted' to free memories; (0.1.0)
	// implement GPU-friendly face array and embed its updates to improve
	// framerate; (0.5.8)
	// move garbage_collection here since we have a better free memory re-allocation
	// strategy; (0.5.9)
	garbage_collection(true, true, true);

	return need_further_refined_;
}

ProgressiveMesh::VertexHandle ViewDependentProgressiveMesh::new_vertex()
{
	if (free_vertex_index_ == -1)
	{
		VertexHandle vh = ProgressiveMesh::new_vertex();
		size_t nV = n_vertices();
		vertex_params_.resize(nV);
		vertex_status_.resize(nV);
		vertex_remark_.resize(nV);
		return vh;
	}
	else
	{
		VertexHandle i0 = free_vertices_[free_vertex_index_--];
		vertex_status_[i0].set_deleted(false);
		return i0;
     }
}

ProgressiveMesh::VertexHandle ViewDependentProgressiveMesh::new_vertex(const Point& p)
{
	return ProgressiveMesh::new_vertex(p);
}

ProgressiveMesh::HalfedgeHandle ViewDependentProgressiveMesh::new_edge()
{
	EdgeHandle eidx = -1;
	if (free_edge_index_ == -1)
	{
		eidx = ProgressiveMesh::new_edge() >> 1;
		edge_status_.resize(n_edges());
	}
	else
	{
		eidx = free_edges_[free_edge_index_--];
		edge_status_[eidx].set_deleted(false);
	}
	return eidx<<1;
}

ProgressiveMesh::HalfedgeHandle ViewDependentProgressiveMesh::new_edge(
	VertexHandle start_vertex_handle, VertexHandle end_vertex_handle)
{
	return ProgressiveMesh::new_edge(start_vertex_handle, end_vertex_handle);
}

ProgressiveMesh::FaceHandle ViewDependentProgressiveMesh::new_face()
{
	FaceHandle f1 = -1;
	if (free_face_index_ == -1)
		f1 = ProgressiveMesh::new_face();
	else
	{
		f1 = free_faces_[free_face_index_--];
		face_status_[f1].set_deleted(false);
	}
	return f1;
}

ProgressiveMesh::FaceHandle ViewDependentProgressiveMesh::new_face(
	VertexHandle i0, VertexHandle i1, VertexHandle i2)
{
	return ProgressiveMesh::new_face(i0, i1, i2);
}

bool ViewDependentProgressiveMesh::force_vsplit(std::list<VParam>::iterator vpit)
{
    const VsData& vs = loadVs(vpit->pos);

	VertexHandle vl = -1, vr = -1;
	HalfedgeHandle v1vl = -1, vrv1 = -1;
	get_active_cuts(vpit, vs, vl, vr, v1vl, vrv1);

	while (vl == vr)
	{
		if (!force_vsplit (vertex_params_[vl]))
            return false;
		get_active_cuts(vpit, vs, vl, vr, v1vl, vrv1);
	}

	return vsplit(vpit, vs, vl, vr, v1vl, vrv1);
}

void ViewDependentProgressiveMesh::ecol(
	std::list<VParam>::iterator vpit, HalfedgeHandle v0v1)
{
	VertexHandle v1 = vpit->rchild->vh;

	collapse(v0v1);
	get_vertex_normal(v1) = vpit->normal;
	vertex_params_[v1] = vpit;
	vpit->vh = v1;

	// mark the free vparams from edge-collapses; (0.5.9)
	free_vplist_[++free_param_index_] = vpit->lchild;
	vpit->lchild = vplist_.end();
	free_vplist_[++free_param_index_] = vpit->rchild;
	vpit->rchild = vplist_.end();
}

void ViewDependentProgressiveMesh::unload_mesh()
{
    kappa_square_ = 0.0;
    tree_id_bits_ = 0;
    base_offset_ = 0;
    data_size_ = 0;
    free_vertex_index_ = -1;
    free_face_index_ = -1;
    free_edge_index_ = -1;
    free_param_index_ = -1;
    is_mesh_loaded_ = false;
    need_refined_ = false;

    free_vertices_.clear();
    free_vertices_.reserve(MAX_ECVS_PER_REFINEMENT);
    
    free_faces_.clear();
    free_faces_.reserve(2*MAX_ECVS_PER_REFINEMENT);

    free_edges_.clear();
    free_edges_.reserve(3*MAX_ECVS_PER_REFINEMENT);

    free_vplist_.clear();
    free_vplist_.reserve(2*MAX_ECVS_PER_REFINEMENT);

    ProgressiveMesh::unload_mesh();

    deactivate();

	is_mesh_loaded_ = false;
}

void ViewDependentProgressiveMesh::calc_tree_id_bits(unsigned int num_roots)
{
	tree_id_bits_ = 0;
	while (num_roots > ((unsigned int) 0x00000001 << tree_id_bits_))
		++tree_id_bits_;

	tree_id_bitshift_ = 32 - tree_id_bits_;
	tree_id_mask_     = 0xFFFFFFFF >> tree_id_bits_;
}

bool ViewDependentProgressiveMesh::qrefine(const VParam& vp)
{
    float tolerance_square = viewing_parameters_.tolerance_square() / 4025.0f;
	float tan_value = tanf (45.0f / 2.0f);
	kappa_square_   = 4.0f * tan_value * tan_value * tolerance_square;

	OpenMesh::Vec3f& p = get_point(vp.vh);
	if (viewing_parameters_.outside_view_frustum(p, vp.radius) == true)
		return false;

	OpenMesh::Vec3f eye_dir = p - viewing_parameters_.eye_pos();
	float  distance         = eye_dir.length();
	float  distance2        = distance * distance;
	float  product_value    = dot(eye_dir, vp.normal);
	if (oriented_away(vp.sin_square, distance2, product_value) == true)
		return false;
	if (screen_space_error(vp.mue_square, vp.sigma_square, distance2, product_value)
		== true)
	{
		return false;
	}
	return true;
}

bool ViewDependentProgressiveMesh::oriented_away(
	float sin_square, float distance_square, float product_value) const
{
	if (product_value > 0 &&
		product_value * product_value > distance_square * sin_square)
	{
		return true;
	}
	return false;
}

bool ViewDependentProgressiveMesh::screen_space_error(
	float mue_square, float sigma_square,
	float distance_square, float product_value) const
{
	if ((mue_square >= kappa_square_ * distance_square) ||
		(sigma_square*(distance_square-product_value*product_value) >=
		 kappa_square_*distance_square*distance_square))
	{
		return false;
	}
	return true;
}

void ViewDependentProgressiveMesh::get_active_cuts(
	std::list<VParam>::iterator vpit, const VsData& vs,
	VertexHandle& vl, VertexHandle& vr,
	HalfedgeHandle& v1vl, HalfedgeHandle& vrv1) const
{
	// Note! We remove the dependencies on VHierarchyNodeIndex to save some CPU
	// time; (0.5.5)
	VertexHandle vv_it = -1;
	vl = -1;
	vr = -1;
	v1vl = -1;
	vrv1 = -1;

	VertexID neighbor;
	VertexID fund_lcut = vs.fund_lcut_index;
	VertexID fund_rcut = vs.fund_rcut_index;

	unsigned int fund_lcut_tree_id = fund_lcut >> (tree_id_bitshift_);
	unsigned int fund_rcut_tree_id = fund_rcut >> (tree_id_bitshift_);
	unsigned int fund_lcut_node_id = fund_lcut &  (tree_id_mask_);
	unsigned int fund_rcut_node_id = fund_rcut &  (tree_id_mask_);

	HalfedgeHandle sheh = get_vertex(vpit->vh).halfedge_handle_;
	HalfedgeHandle  heh = sheh;
	do
	{
		const Halfedge& hehe = get_halfedge(heh);
		vv_it = hehe.vertex_handle_;
		neighbor = vertex_params_[vv_it]->vid;

		unsigned int  neighbor_tree_id = neighbor >> (tree_id_bitshift_);
		unsigned int  neighbor_node_id = neighbor &  (tree_id_mask_);

		// finds the halfedges (v1, vl), (vr, v1) here to save one neighborhood
		// traverse; (0.5.5)

		// Find the left cut vertex: neighbor_node_id is an ancestor of the
		// fund_lcut_node_id (and they belong to the same tree) - this is necessary
		// because the actual left cut vertex might not be active. However, a
		// neighbor can still be considered a valid cut if it is an ancestor of the
		// actual left cut node.
		if (!is_valid_vertex(vl)
			&& neighbor_tree_id == fund_lcut_tree_id
			&& neighbor_node_id <= fund_lcut_node_id)
		{
			unsigned int lcid = fund_lcut_node_id;
			while (lcid > neighbor_node_id)
				lcid >>= 1;
			if (neighbor_node_id == lcid) // v1vl = find_halfedge(v1, vl);
				vl = vv_it; v1vl = heh;
		}

		if (is_valid_vertex(vl) && is_valid_vertex(vr))
			break;

		if (!is_valid_vertex(vr)
			&& neighbor_tree_id == fund_rcut_tree_id
			&& neighbor_node_id <= fund_rcut_node_id)
		{
			unsigned int rcid = fund_rcut_node_id;
			while (rcid > neighbor_node_id)
				rcid >>= 1;
			if (neighbor_node_id == rcid) // vrv1 = find_halfedge(vr, v1);
				vr = vv_it; vrv1 = get_opposite_halfedge_handle(heh);
		}

		if (is_valid_vertex(vl) && is_valid_vertex(vr))
			break;

		heh = get_opposite_halfedge(heh).next_halfedge_handle_;
	} while (heh != sheh);
}

bool ViewDependentProgressiveMesh::vsplit(
	std::list<VParam>::iterator vpit, const VsData& vs,
	VertexHandle vl, VertexHandle vr,
	HalfedgeHandle v1vl, HalfedgeHandle vrv1)
{
    if (n_vertices() >= Constants::MAX_VERTICES)
        return false;

	VHierarchyNodeIndex   node_index = VHierarchyNodeIndex(vpit->vid);
	VHierarchyNodeIndex lchild_index =
		VHierarchyNodeIndex(
			node_index.tree_id(tree_id_bits_),
			2*node_index.node_id(tree_id_bits_),
			tree_id_bits_);
	VHierarchyNodeIndex rchild_index =
		VHierarchyNodeIndex(
			node_index.tree_id(tree_id_bits_),
			2*node_index.node_id(tree_id_bits_)+1,
			tree_id_bits_);

	VParam lvp, rvp;

	lvp.vh            = new_vertex(vs.p);
	lvp.pos           = vs.l_pos;
	lvp.vid           = lchild_index.value();
	lvp.normal        = vs.l_normal;
	lvp.radius        = vs.l_radius;
	lvp.sin_square    = vs.l_sin_square;
	lvp.mue_square    = vs.l_mue_square;
	lvp.sigma_square  = vs.l_sigma_square;
	lvp.parent        = vpit;
	lvp.lchild        = vplist_.end();
	lvp.rchild        = vplist_.end();

	rvp.vh            = vpit->vh;
	rvp.pos           = vs.r_pos;
	rvp.vid           = rchild_index.value();
	rvp.normal        = vs.r_normal;
	rvp.radius        = vs.r_radius;
	rvp.sin_square    = vs.r_sin_square;
	rvp.mue_square    = vs.r_mue_square;
	rvp.sigma_square  = vs.r_sigma_square;
	rvp.parent        = vpit;
	rvp.lchild        = vplist_.end();
	rvp.rchild        = vplist_.end();

	// re-allocate the free vparams to vertex-splits; (0.5.9)
	if (free_param_index_ == -1)
		vpit->lchild = vplist_.insert(vplist_.end(), lvp);
	else
	{
		vpit->lchild = free_vplist_[free_param_index_--];
		(*(vpit->lchild)) = lvp;
	}

	if (free_param_index_ == -1)
		vpit->rchild = vplist_.insert(vplist_.end(), rvp);
	else
	{
		vpit->rchild = free_vplist_[free_param_index_--];
		(*(vpit->rchild)) = rvp;
	}

	VertexHandle v0 = lvp.vh, v1 = rvp.vh;
	ProgressiveMesh::vertex_split(v0, v1, vl, vr, v1vl, vrv1);
	get_vertex_normal(v0) = lvp.normal;
	get_vertex_normal(v1) = rvp.normal;
	vertex_params_[v0] = vpit->lchild;
	vertex_params_[v1] = vpit->rchild;
    
    return true;
}

bool ViewDependentProgressiveMesh::is_collapse_ok(
	VertexHandle v0, VertexHandle v1, HalfedgeHandle& v0v1)
{
    // are vertices already deleted ?
    if (vertex_status_[v0].deleted() || vertex_status_[v1].deleted())
        return false;

    // Note! We find the halfedge (v0, v1) here and remark v0's
    // neighbors before intersection tests; (0.5.6)
    HalfedgeHandle sheh = get_vertex(v0).halfedge_handle_;
    HalfedgeHandle  heh = sheh;
    VertexHandle     vh = -1;
    do {
        vh = get_halfedge(heh).vertex_handle_;
        vertex_remark_[vh] = v0;
        if (vh == v1)
            v0v1 = heh;
        heh = get_opposite_halfedge(heh).next_halfedge_handle_;
    } while (heh != sheh);

    Halfedge&   he = get_halfedge(v0v1); // v0v1; v1 = he.vertex_handle_;
    Halfedge&  ohe = get_opposite_halfedge(v0v1); // v1v0; v0 = ohe.vertex_handle_;

    VertexHandle vl = -1, vr = -1;

    // the edges v1-vl and vl-v0 must not be both boundary edges
    if (is_valid_face(he.face_handle_))
    {
        HalfedgeHandle ph = he.prev_halfedge_handle_;
        HalfedgeHandle nh = he.next_halfedge_handle_;
        Halfedge&    ophe = get_opposite_halfedge(ph);
        Halfedge&    onhe = get_opposite_halfedge(nh);
        if (!is_valid_face(ophe.face_handle_) && 
            !is_valid_face(onhe.face_handle_)) return false;
        vl = ophe.vertex_handle_;
    }

    // the edges v0-vr and vr-v1 must not be both boundary edges
    if (is_valid_face(ohe.face_handle_))
    {
        HalfedgeHandle ph = ohe.prev_halfedge_handle_;
        HalfedgeHandle nh = ohe.next_halfedge_handle_;
        Halfedge&    ophe = get_opposite_halfedge(ph);
        Halfedge&    onhe = get_opposite_halfedge(nh);
        if (!is_valid_face(ophe.face_handle_) && 
            !is_valid_face(onhe.face_handle_)) return false;
        vr = ophe.vertex_handle_;
    }
    
    // if vl and vr are equal or both invalid -> fail
    if (vl == vr) return false;
    
    // edge between two boundary vertices should be a boundary edge
    if ( is_boundary_vertex(v0) && is_boundary_vertex(v1) &&
        is_valid_face(he.face_handle_) && is_valid_face(ohe.face_handle_))
        return false;

    // Note! We standardize the namings and routines and optimize a little; (0.5.6)
    sheh = get_vertex(v1).halfedge_handle_;
    heh = sheh;
    do {
        vh = get_halfedge(heh).vertex_handle_;
        heh = get_opposite_halfedge(heh).next_halfedge_handle_;
        if (vertex_remark_[vh] == v0 && vh != vl && vh != vr)
            return false;
    } while (heh != sheh);

    return true; // passed all tests
}

void ViewDependentProgressiveMesh::collapse(HalfedgeHandle hh)
{
    HalfedgeHandle  v0v1 = hh;
    Halfedge&  he01 = get_halfedge(v0v1);
    HalfedgeHandle  v1v0 = get_opposite_halfedge_handle(v0v1);
    Halfedge&  he10 = get_halfedge(v1v0);
    
    HalfedgeHandle  hn = he01.next_halfedge_handle_;
    Halfedge&      hne = get_halfedge(hn);
    HalfedgeHandle  hp = he01.prev_halfedge_handle_;
    Halfedge&      hpe = get_halfedge(hp);
    HalfedgeHandle  on = he10.next_halfedge_handle_;
    Halfedge&      one = get_halfedge(on);
    HalfedgeHandle  op = he10.prev_halfedge_handle_;
    Halfedge&      ope = get_halfedge(op);
    
    FaceHandle   fh = he01.face_handle_;
    FaceHandle   fo = he10.face_handle_;
    VertexHandle v1 = he01.vertex_handle_;
    VertexHandle v0 = he10.vertex_handle_;
    
    // Note! We do not always need to adjust the outgoing (boundary) halfedge of v1;
    Vertex& vv0 = get_vertex(v0);
    Vertex& vv1 = get_vertex(v1); 
    HalfedgeHandle heh = v0v1;
    FaceHandle df = -1;
    FaceHandle pdf = he01.face_handle_;
    do {
        Halfedge&  hehe = get_halfedge(heh);
        Halfedge& ohehe = get_opposite_halfedge(heh);
        ohehe.vertex_handle_ = v1;

        df = ohehe.face_handle_; int vd = hehe.vertex_handle_;
        if ( df >= 0) draw_faces_[ df*3  ] = v1; // ohehe.vertex_handle_;
        if (pdf >= 0) draw_faces_[pdf*3+1] = vd; //  hehe.vertex_handle_;
        if ( df >= 0) draw_faces_[ df*3+2] = vd; //  hehe.vertex_handle_;
        pdf = df;

        heh = ohehe.next_halfedge_handle_;
    } while (heh != v0v1);

    // adjust outgoing halfedge of v1
    Halfedge& vhe0 = get_halfedge(vv0.halfedge_handle_);
    if (vv0.halfedge_handle_ != v0v1 && vhe0.face_handle_ == -1)
        vv1.halfedge_handle_ = vv0.halfedge_handle_; 
    if (vv1.halfedge_handle_ == v1v0)
        vv1.halfedge_handle_ = hn;

    hpe.next_halfedge_handle_ = hn;
    hne.prev_halfedge_handle_ = hp;
    ope.next_halfedge_handle_ = on;
    one.prev_halfedge_handle_ = op;

    if (is_valid_face(fh))
        get_face(fh).halfedge_handle_ = hn;
    if (is_valid_face(fo))
        get_face(fo).halfedge_handle_ = on;

    // Note! We mark the free edge and vertex from edge-collapses; (0.5.9)
    get_vertex(v0).halfedge_handle_ = -1; // invalidate(), set_isolated(v0);
    edge_status_[v0v1>>1].set_deleted(true);
    free_edges_[++free_edge_index_] = (v0v1>>1);
    vertex_status_[v0].set_deleted(true);
    free_vertices_[++free_vertex_index_] = v0;
    vertex_params_[v0] = vplist_.end();
    vertex_remark_[v0] = -1;

    // remove loops, collapse_loop(next_halfedge_handle(hn)) 
    if (hne.next_halfedge_handle_ == hp && hpe.prev_halfedge_handle_ == hn)
        collapse_loop(hp, hpe, hn, hne);

    // remove loops, collapse_loop(on)
    if (one.next_halfedge_handle_ == op && ope.prev_halfedge_handle_ == on)
        collapse_loop(on, one, op, ope);
}

void ViewDependentProgressiveMesh::collapse_loop(
	HalfedgeHandle hp, Halfedge& hpe,
	HalfedgeHandle hn, Halfedge& hne)
{
	HalfedgeHandle h0 = hp;
	Halfedge&     he0 = hpe;
	HalfedgeHandle h1 = hn;
	Halfedge&     he1 = hne;

	HalfedgeHandle  o0 = get_opposite_halfedge_handle(h0);
	Halfedge&      oe0 = get_halfedge(o0);
	HalfedgeHandle op0 = oe0.prev_halfedge_handle_;
	Halfedge&     ope0 = get_halfedge(op0);
	HalfedgeHandle on0 = oe0.next_halfedge_handle_;
	Halfedge&     one0 = get_halfedge(on0);

	VertexHandle v0 = he0.vertex_handle_;
	VertexHandle v1 = he1.vertex_handle_;
	FaceHandle   fh = he0.face_handle_;
	FaceHandle   fo = oe0.face_handle_;

	he1.next_halfedge_handle_ = on0;
	one0.prev_halfedge_handle_ =  h1;
	he1.prev_halfedge_handle_ = op0;
	ope0.next_halfedge_handle_ =  h1;
	he1.face_handle_ = fo;

	if (get_vertex(v0).halfedge_handle_ == o0)
		get_vertex(v0).halfedge_handle_ = h1;
	if (get_vertex(v1).halfedge_handle_ == h0)
		get_vertex(v1).halfedge_handle_ = on0;

	if (is_valid_face(fo) && get_face(fo).halfedge_handle_ == o0)
		get_face(fo).halfedge_handle_ = h1;

	// mark the free face and edge from edge-collapses; (0.5.9)
	if (is_valid_face(fh))
	{
		get_face(fh).halfedge_handle_ = -1;
		face_status_[fh].set_deleted(true);
		free_faces_[++free_face_index_] = fh;
	}

	edge_status_[edge_handle(h0)].set_deleted(true);
	free_edges_[++free_edge_index_] = edge_handle(h0);
}

static bool idxcomp (int i, int j) { return (i>j); }

void ViewDependentProgressiveMesh::garbage_collection(bool v, bool e, bool f)
{
    int i0, i1, nV(n_vertices()), nE(n_edges()), nF(n_faces());
    
    // remove deleted vertices
    if (v && n_vertices() > 0)
    {
        // sort the free vertex pointers before traverse them; (0.5.9)
        sort(free_vertices_.begin(),
             free_vertices_.begin()+free_vertex_index_+1,
             &MobiMesh::idxcomp);
        i0 = 0;
        i1 = nV-1;
        while (1)
        {
            // find 1st deleted and last un-deleted, and then swap them; (0.1.0)
            // mark the free vertices from edge-collapses to speed up garbage_collection; (0.5.9)
            if (free_vertex_index_ == -1)
                break;
            i0 = free_vertices_[free_vertex_index_--];

            while ( vertex_status_[i1].deleted() && i0 < i1 )
                --i1;

            if (i0 >= i1)
                break;

            // update handles of halfedges
            HalfedgeHandle sheh = get_vertex(i1).halfedge_handle_;
            Halfedge& shehe = get_halfedge(sheh);
            HalfedgeHandle heh = sheh;
            FaceHandle df = -1;
            FaceHandle pdf = shehe.face_handle_;
            do
            {
                Halfedge&  hehe = get_halfedge(heh);
                Halfedge& ohehe = get_opposite_halfedge(heh);
                ohehe.vertex_handle_ = i0;

                // embed the updates of GPU-friendly 'draw_faces_' into garbage_collection; (0.5.8)
                // fix the segmentation fault while updating draw_faces_ on neptune_nonmanifold.vpm; (0.5.9)
                df = ohehe.face_handle_;
                VertexHandle vd = hehe.vertex_handle_;
                if ( df >= 0) draw_faces_[ df*3  ] = i0; // ohehe.vertex_handle_;
                if (pdf >= 0) draw_faces_[pdf*3+1] = vd; //  hehe.vertex_handle_;
                if ( df >= 0) draw_faces_[ df*3+2] = vd; //  hehe.vertex_handle_;
                pdf = df;

                heh = ohehe.next_halfedge_handle_;
            } while (heh != sheh);

            std::swap(     vertices_[i0],      vertices_[i1]);
            std::swap(       points_[i0],        points_[i1]);
            std::swap(vertex_normal_[i0], vertex_normal_[i1]);
            std::swap(vertex_params_[i0], vertex_params_[i1]);
            std::swap(vertex_status_[i0], vertex_status_[i1]);
            std::swap(vertex_remark_[i0], vertex_remark_[i1]);

            // update VertexHandle links since 'garbage_collection()' will shuffle them;
            vertex_params_[i0]->vh = i0;
        };
        nV = vertex_status_[i1].deleted() ? i1 : i1+1;
        vertices_.resize(nV);
        points_.resize(nV);
        vertex_normal_.resize(nV);
        vertex_params_.resize(nV);
        vertex_status_.resize(nV);
        vertex_remark_.resize(nV);
    }
    // remove deleted faces
    if (f && n_faces() > 0)
    {
        // Note! We need to sort the free face pointers before traverse them; (0.5.9)
        sort(free_faces_.begin(),
             free_faces_.begin()+free_face_index_+1,
             &MobiMesh::idxcomp);
        i0 = 0;
        i1 = nF-1;
        while (1)
        {
            // find 1st deleted and last un-deleted, and then swap them; (0.1.0)
            // mark the free faces from edge-collapses to speed up garbage_collection; (0.5.9)
            if (free_face_index_ == -1)
                break;
            i0 = free_faces_[free_face_index_--];

            while ( face_status_[i1].deleted() && i0 < i1 )
                --i1;

            if (i0 >= i1)
                break;

            // update handles of halfedges
            HalfedgeHandle fh = get_face(i1).halfedge_handle_;
            Halfedge&  fhe = get_halfedge(fh);
            fhe.face_handle_ = i0;

            HalfedgeHandle pfh =   fhe.prev_halfedge_handle_;
            Halfedge& pfhe = edges_[pfh >> 1].halfedges_[pfh & 1];
            pfhe.face_handle_ = i0;

            HalfedgeHandle nfh =   fhe.next_halfedge_handle_;
            Halfedge& nfhe = edges_[nfh >> 1].halfedges_[nfh & 1];
            nfhe.face_handle_ = i0;
            
            std::swap(      faces_[i0],       faces_[i1]);
            std::swap(face_normal_[i0], face_normal_[i1]);
            std::swap(face_status_[i0], face_status_[i1]);

            // swap the elements of GPU-friendly 'draw_faces_' in garbage_collection; (0.5.8)
            draw_faces_[i0*3  ] = pfhe.vertex_handle_;
            draw_faces_[i0*3+1] =  fhe.vertex_handle_;
            draw_faces_[i0*3+2] = nfhe.vertex_handle_;
        };
        nF = face_status_[i1].deleted() ? i1 : i1+1;
        faces_.resize(nF);
        face_normal_.resize(nF);
        face_status_.resize(nF);
        
        draw_faces_.resize(nF*3);
    }
    // remove deleted edges
    if (e && n_edges() > 0)
    {
        // Note! We need to sort the free edge pointers before traverse them; (0.5.9)
        sort(free_edges_.begin(),
             free_edges_.begin()+free_edge_index_+1,
             &MobiMesh::idxcomp);
        i0 = 0;
        i1 = nE-1;
        while (1)
        {
            // find 1st deleted and last un-deleted, and then swap them; (0.1.0)
            // mark the free edges from edge-collapses to speed up garbage_collection; (0.5.9)
            if (free_edge_index_ == -1)
                break;
            i0 = free_edges_[free_edge_index_--];         

            while ( edge_status_[i1].deleted() && i0 < i1 )
                --i1;

            if (i0 >= i1)
                break;
            
            HalfedgeHandle h0 = 2*i0;
            HalfedgeHandle h1 = 2*i0+1;
            Halfedge&  he0 = edges_[i1].halfedges_[0];
            Halfedge&  he1 = edges_[i1].halfedges_[1];

            // update handles of vertices
            Vertex& v0 = get_vertex(he0.vertex_handle_);
            Vertex& v1 = get_vertex(he1.vertex_handle_);
            if (v0.halfedge_handle_ == 2*i1+1)
                v0.halfedge_handle_ = h1;
            if (v1.halfedge_handle_ == 2*i1  )
                v1.halfedge_handle_ = h0;

            // update handles of halfedges (next/prev)
            HalfedgeHandle ph0 = he0.prev_halfedge_handle_;
            HalfedgeHandle nh0 = he0.next_halfedge_handle_;
            get_halfedge(ph0).next_halfedge_handle_ = h0;
            get_halfedge(nh0).prev_halfedge_handle_ = h0;
            HalfedgeHandle ph1 = he1.prev_halfedge_handle_;
            HalfedgeHandle nh1 = he1.next_halfedge_handle_;
            get_halfedge(ph1).next_halfedge_handle_ = h1;
            get_halfedge(nh1).prev_halfedge_handle_ = h1;

            // update handles of faces
            Face&   f0 = get_face(he0.face_handle_);
            Face&   f1 = get_face(he1.face_handle_);
            if (f0.halfedge_handle_ == 2*i1  )
                f0.halfedge_handle_ = h0;
            if (f1.halfedge_handle_ == 2*i1+1)
                f1.halfedge_handle_ = h1;

            std::swap(      edges_[i0],       edges_[i1]);
            std::swap(edge_status_[i0], edge_status_[i1]);
        };
        nE = edge_status_[i1].deleted() ? i1 : i1+1;
        edges_.resize(nE);
        edge_status_.resize(nE);
    }

    // Note! We reset all the free memory pointers here just for sure; (0.5.9)
    free_vertex_index_ = -1;
    free_edge_index_ = -1;
    free_face_index_ = -1;

    // Note! We need to do garbage_collection on free_vplist_; (0.6.2)
    while (free_param_index_ != -1 )
        vplist_.erase(free_vplist_[free_param_index_--]);
}

inline ViewDependentProgressiveMesh::VsData
ViewDependentProgressiveMesh::loadVs(unsigned int data_index)
{
    if (mmapped_)
        return (VsData &) * (data_ + base_offset_ + (data_index-1)*data_size_);
    else
    {
        VsData data;
        pread(fd_, &data, sizeof(VsData), (base_offset_ + (data_index-1)*data_size_));
        return data;
    }
}

bool ViewDependentProgressiveMesh::is_ancestor_VertexID(
	VertexID ancestor, VertexID descendent) const
{
	unsigned int   ancestor_tree_id =   ancestor >> (32 - tree_id_bits_);
	unsigned int descendent_tree_id = descendent >> (32 - tree_id_bits_);

	if (ancestor_tree_id != descendent_tree_id)
		return false;

	unsigned int ancestor_node_id =
		ancestor & ((unsigned int) 0xFFFFFFFF >> tree_id_bits_);
	unsigned int descendent_node_id =
		descendent & ((unsigned int) 0xFFFFFFFF >> tree_id_bits_);

	if (ancestor_node_id  > descendent_node_id)
		return false;

	while (descendent_node_id > 0)
	{
		if (ancestor_node_id == descendent_node_id)
			return true;
		descendent_node_id >>= 1;
	}
	return false;
}

}
