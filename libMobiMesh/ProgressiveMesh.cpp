/*
 *  ProgressiveMesh.cpp
 *  MobiMesh
 *
 *  Created by Ong Yuh Shin on 3/9/11.
 *  Copyright 2011 National University of Singapore. All rights reserved.
 *
 */

#include <fstream>
#include <sys/types.h>

#include "MobiMeshConstants.h"
#include "FileReadErrorException.h"
#include "ProgressiveMesh.h"
#include "InvalidFileFormatException.h"

namespace MobiMesh {

///////////// Constructors

ProgressiveMesh::ProgressiveMesh(const std::string& filename)
{
	load_mesh(filename);
}


///////////// Loading and saving mesh data

void ProgressiveMesh::load_mesh(const std::string& filename)
{
	char            fileformat[10];
	OpenMesh::Vec3f p;
	unsigned int    v1,vl,vr,fvi[3];

	std::ifstream ifs(filename.c_str(), std::ios::binary);

	if (!ifs)
		throw FileReadErrorException("ifstream read error");
	
	// first 8 bytes/chars should be file type "ProgMesh"
	ifs.read(fileformat, 8); fileformat[8]='\0';
	if (std::string(fileformat) != std::string("ProgMesh"))
	{
		ifs.close();
		throw InvalidFileFormatException(
			"invalid file format", // error msg
			"ProgMesh",            // expected format
			fileformat);           // received format
	}

	// read number of vertices & faces in base mesh and number of vertex splits
	ifs.read( (char*)&n_base_vertices_, sizeof(n_base_vertices_) );
	ifs.read( (char*)&n_base_faces_   , sizeof(n_base_faces_) );
	ifs.read( (char*)&n_details_      , sizeof(n_details_) );

	// read the positions of the vertices of the base mesh
	for (size_t i=0; i<n_base_vertices_; ++i)
	{
		ifs.read( (char*)&p , sizeof(p) );
		new_vertex(p);
	}

	// read the faces (location of its vertices)
	for (size_t i=0; i<n_base_faces_; ++i)
	{
		ifs.read( (char*)&fvi[0] , sizeof(fvi[0]) );
		ifs.read( (char*)&fvi[1] , sizeof(fvi[1]) );
		ifs.read( (char*)&fvi[2] , sizeof(fvi[2]) );
		new_face(fvi[0],fvi[1],fvi[2]);
	}

	// read in subsequent vertex split data and perform the split until we have
	// too many vertices (> n_max_vertices_)
	unsigned int n_max_vertices_ = n_base_vertices_ + n_details_;
	n_max_vertices_ = (n_max_vertices_ > Constants::MAX_VERTICES) ?
		Constants::MAX_VERTICES : n_max_vertices_;
	for (size_t i=n_base_vertices_; i<n_max_vertices_; ++i)
	{
		ifs.read( (char*)&p  , sizeof(p)  );
		ifs.read( (char*)&v1 , sizeof(v1) );
		ifs.read( (char*)&vl , sizeof(vl) );
		ifs.read( (char*)&vr , sizeof(vr) );

		// Note! We need to meet the requirements of find_halfedge; (0.6.2)
		HalfedgeHandle v1vl = (vl == -1) ? -1 : find_halfedge(v1, vl);
		HalfedgeHandle vrv1 = (vr == -1) ? -1 : find_halfedge(vr, v1);
		vertex_split(new_vertex(p), v1, vl, vr, v1vl, vrv1);
	}

	ifs.close();
	
	update_normals();
}
    
void ProgressiveMesh::unload_mesh()
{
    draw_faces_.clear();

	points_.clear();
	vertex_normal_.clear();

    vertices_.clear();
	edges_.clear();
	faces_.clear();
}

ProgressiveMesh::VertexHandle ProgressiveMesh::new_vertex()
{
	vertices_.push_back(Vertex());
	points_.push_back(IMesh::Point(0,0,0));
	size_t num_vertices = n_vertices();
	vertex_normal_.resize(num_vertices);
	return num_vertices-1;
}

ProgressiveMesh::VertexHandle ProgressiveMesh::new_vertex(
	const IMesh::Point& p)
{
	VertexHandle new_vertex_handle = new_vertex();
	points_[new_vertex_handle] = p;
	return new_vertex_handle;
}

ProgressiveMesh::HalfedgeHandle ProgressiveMesh::new_edge()
{
	edges_.push_back(Edge());
	return (n_edges()-1)<<1;
}

ProgressiveMesh::HalfedgeHandle ProgressiveMesh::new_edge(
	VertexHandle start_vertex_handle, VertexHandle end_vertex_handle)
{
	HalfedgeHandle heh = new_edge();
	get_halfedge(heh).vertex_handle_          =   end_vertex_handle;
	get_opposite_halfedge(heh).vertex_handle_ = start_vertex_handle;
	return heh;
}

ProgressiveMesh::FaceHandle ProgressiveMesh::new_face()
{
	faces_.push_back(Face());
	draw_faces_.resize(n_faces()*3);
	return n_faces()-1;
}

ProgressiveMesh::FaceHandle ProgressiveMesh::new_face(
	VertexHandle i0, VertexHandle i1, VertexHandle i2)
{
	VertexHandle vertex_handles[3] = { i0, i1, i2 };
	uint vhs_size = 3;

	VertexHandle        vh;
	uint                i, ii, n(vhs_size), vid;
	std::vector<int>    halfedge_handles(n);
	std::vector<bool>   is_new(n), needs_adjust(n, false);
	HalfedgeHandle      inner_next, inner_prev,
	                    outer_next, outer_prev,
	                    boundary_next, boundary_prev,
	                    patch_start, patch_end;

	// cache for set_next_halfedge and vertex' set_halfedge
	typedef NextCacheEntryT<int>           NextCacheEntry;
	typedef std::vector<NextCacheEntry>    NextCache;
	NextCache next_cache;
	next_cache.reserve(3*n);
	assert (n > 2); // don't allow degenerated faces

	// test for topological errors
	// ensure that newly added face does not break manifold surface of
	// the mesh, meaning we cannot add faces that contain non-boundary
	// vertices or non-boundary edges.
	for (i=0, ii=1; i<n; ++i, ++ii, ii%=n)
	{
		if ( !is_boundary_vertex(vertex_handles[i]) )
		{
			std::cerr << "ProgressiveMesh::new_face: complex vertex\n";
			return -1;
		}
		halfedge_handles[i] =
			find_halfedge(vertex_handles[i], vertex_handles[ii]);
		is_new[i] = !is_valid_halfedge(halfedge_handles[i]);
		if (!is_new[i] && !is_boundary_halfedge(halfedge_handles[i]))
		{
			std::cerr << "ProgressiveMesh::new_face: complex edge\n";
			return -1;
		}
	}

	// re-link patches if necessary
	for (i=0, ii=1; i<n; ++i, ++ii, ii%=n)
	{
		if (!is_new[i] && !is_new[ii])
		{
			inner_prev = halfedge_handles[i];
			inner_next = halfedge_handles[ii];

			if (next_halfedge_handle(inner_prev) != inner_next)
			{
				// here comes the ugly part... we have to relink a whole patch
				// search a free gap
				// free gap will be between boundary_prev and boundary_next
				outer_prev = get_opposite_halfedge_handle(inner_next);
				outer_next = get_opposite_halfedge_handle(inner_prev);
				boundary_prev = outer_prev;
				do
				{
					boundary_prev = get_opposite_halfedge_handle(
						next_halfedge_handle(boundary_prev));
				} while (!is_boundary_halfedge(boundary_prev) ||
						boundary_prev==inner_prev);
				boundary_next = next_halfedge_handle(boundary_prev);
				assert(is_boundary_halfedge(boundary_prev));
				assert(is_boundary_halfedge(boundary_next));
				// ok ?
				if (boundary_next == inner_next)
				{
					std::cerr << "ProgressiveMesh::new_face: patch re-linking failed\n";
					return -1;
				}
				// other halfedges' handles
				patch_start = next_halfedge_handle(inner_prev);
				patch_end   = prev_halfedge_handle(inner_next);
				// relink
				next_cache.push_back(NextCacheEntry(boundary_prev, patch_start));
				next_cache.push_back(NextCacheEntry(patch_end, boundary_next));
				next_cache.push_back(NextCacheEntry(inner_prev, inner_next));
			}
		}
	}

	// create missing edges
	for (i=0, ii=1; i<n; ++i, ++ii, ii%=n)
		if (is_new[i])
			halfedge_handles[i] = new_edge(vertex_handles[i], vertex_handles[ii]);

	// create the face
	FaceHandle fh = new_face();
	set_halfedge_handle_fh(fh, halfedge_handles[n-1]);

	// setup halfedges
	for (i=0, ii=1; i<n; ++i, ++ii, ii%=n)
	{
		vh         =  vertex_handles[ii];
		inner_prev = halfedge_handles[i];
		inner_next = halfedge_handles[ii];

		vid = 0;
		if (is_new[i])  vid |= 1;
		if (is_new[ii]) vid |= 2;

		if (vid)
		{
			outer_prev = get_opposite_halfedge_handle(inner_next);
			outer_next = get_opposite_halfedge_handle(inner_prev);

			// set outer links
			switch (vid)
			{
				case 1: // prev is new, next is old
					boundary_prev = prev_halfedge_handle(inner_next);
					next_cache.push_back(NextCacheEntry(boundary_prev, outer_next));
					set_halfedge_handle_vh(vh, outer_next);
					break;
				case 2: // next is new, prev is old
					boundary_next = next_halfedge_handle(inner_prev);
					next_cache.push_back(NextCacheEntry(outer_prev, boundary_next));
					set_halfedge_handle_vh(vh, boundary_next);
					break;
				case 3: // both are new
					if (halfedge_handle_vh(vh) == -1)
					{
						set_halfedge_handle_vh(vh, outer_next);
						next_cache.push_back(NextCacheEntry(outer_prev, outer_next));
					}
					else
					{
						boundary_next = halfedge_handle_vh(vh);
						boundary_prev = prev_halfedge_handle(boundary_next);
						next_cache.push_back(NextCacheEntry(boundary_prev, outer_next));
						next_cache.push_back(NextCacheEntry(outer_prev, boundary_next));
					}
					break;
			}
			// set inner link
			next_cache.push_back(NextCacheEntry(inner_prev, inner_next));
		}
		else
			needs_adjust[ii] = (halfedge_handle_vh(vh) == inner_next);
		
		// set face handle
		set_face_handle(halfedge_handles[i], fh);
	}

	// process next halfedge cache
	NextCache::const_iterator ncIt(next_cache.begin()), ncEnd(next_cache.end());
	for (; ncIt != ncEnd; ++ncIt)
		set_next_halfedge_handle(ncIt->first, ncIt->second);

	// adjust vertices' halfedge handle
	for (i=0; i<n; ++i)
		if (needs_adjust[i])
			adjust_outgoing_halfedge(vertex_handles[i]);

	return fh;
}

ProgressiveMesh::HalfedgeHandle ProgressiveMesh::find_halfedge(
	VertexHandle start_vertex_handle,
	VertexHandle end_vertex_handle) const
{
	assert(start_vertex_handle != -1 && end_vertex_handle != -1);
	
	HalfedgeHandle shh = halfedge_handle_vh(start_vertex_handle);
	HalfedgeHandle  hh = shh;
	do
	{
		if (!is_valid_halfedge(hh))
			break;
		if ( get_halfedge(hh).vertex_handle_ == end_vertex_handle )
			return hh;
		hh = get_halfedge(
			get_opposite_halfedge_handle(hh)).next_halfedge_handle_;
	} while (hh != shh);
	
	return -1;
}

ProgressiveMesh::HalfedgeHandle ProgressiveMesh::vertex_split(
	VertexHandle v0, VertexHandle v1, VertexHandle vl, VertexHandle vr,
	HalfedgeHandle v1vl, HalfedgeHandle vrv1)
{
	HalfedgeHandle vlv1 = -1, v1vr = -1, v0v1 = -1;

	// build loop from halfedge vr->v1, insert_loop(vrv1);
	if (vl != -1)
		vlv1 = build_loop(v1, vl, v1vl);
	// build loop from halfedge vr->v1, insert_loop(vrv1);
	if (vr != -1)
		v1vr = build_loop(vr, v1, vrv1);

	// handle boundary cases
	if (vl == -1)
		vlv1 = prev_halfedge_handle(halfedge_handle_vh(v1));
	if (vr == -1)
		vrv1 = prev_halfedge_handle(halfedge_handle_vh(v1));

	// split vertex v1 into edge v0v1, v0v1 = insert_edge(v0, vlv1, vrv1);
	v0v1 = new_edge(v0, v1);
	HalfedgeHandle v1v0 = get_opposite_halfedge_handle(v0v1);

	Halfedge&      he01 = get_halfedge(v0v1);
	Halfedge&      he10 = get_halfedge(v1v0);
	Halfedge&      hel1 = get_halfedge(vlv1);
	Halfedge&      her1 = get_halfedge(vrv1);
	HalfedgeHandle nhl1 = hel1.next_halfedge_handle_;
	HalfedgeHandle nhr1 = her1.next_halfedge_handle_;
	Halfedge&     nhel1 = get_halfedge(nhl1);
	Halfedge&     nher1 = get_halfedge(nhr1);

	he01.next_halfedge_handle_ = nhl1;
	nhel1.prev_halfedge_handle_ = v0v1;
	hel1.next_halfedge_handle_ = v0v1;
	he01.prev_halfedge_handle_ = vlv1;
	he10.next_halfedge_handle_ = nhr1;
	nher1.prev_halfedge_handle_ = v1v0;
	her1.next_halfedge_handle_ = v1v0;
	he10.prev_halfedge_handle_ = vrv1;

	he01.face_handle_ = hel1.face_handle_;
	he10.face_handle_ = her1.face_handle_;

	if (he01.face_handle_ != -1)
		faces_[he01.face_handle_].halfedge_handle_ = v0v1;
	if (he10.face_handle_ != -1)
		faces_[he10.face_handle_].halfedge_handle_ = v1v0;

	// we do not always need to adjust the outgoing (boundary) halfedge of v1.
	Vertex&          vv0 = get_vertex(v0);
	Vertex&          vv1 = get_vertex(v1);
	HalfedgeHandle hehv1 = vv1.halfedge_handle_;
    vv0.halfedge_handle_ = v0v1;
	HalfedgeHandle   heh = v0v1;
    FaceHandle df = -1, pdf = he01.face_handle_;
	do
	{
		Halfedge&  hehe = get_halfedge(heh);
		Halfedge& ohehe = get_opposite_halfedge(heh);
		if (heh == hehv1)
			vv1.halfedge_handle_ = v1v0; // adjust outgoing halfedge of v1
		if (hehe.face_handle_ == -1)
			vv0.halfedge_handle_ = heh;  // adjust outgoing halfedge of v0
		ohehe.vertex_handle_ = v0;

		// embed the updates of GPU-friendly 'draw_faces_' into
		// vertex-splits; (0.5.8)
		// fix the segmentation fault while updating draw_faces_ on
		// neptune_nonmanifold.vpm; (0.5.9)
		df = ohehe.face_handle_;
		VertexHandle vd = hehe.vertex_handle_;
		if ( df >= 0)
			draw_faces_[ df*3  ] = v0; // ohehe.vertex_handle_;
		if (pdf >= 0)
			draw_faces_[pdf*3+1] = vd; //  hehe.vertex_handle_;
		if ( df >= 0)
			draw_faces_[ df*3+2] = vd; //  hehe.vertex_handle_;
		pdf = df;

		heh = ohehe.next_halfedge_handle_;
	} while (heh != v0v1);

	return v0v1;
}

void ProgressiveMesh::update_normals()
{
	int vertices_size = vertices_.size();

	// Keeps track of the number of faces adjacent to each vertex
	std::vector<unsigned short> vertex_faces(vertices_size, 0);

	// Set all vertex normals to 0
	for (int vi = 0; vi < vertices_size; ++vi)
	{
		vertex_normal_[vi] = Normal(0,0,0);
	}

	int faces_size = faces_.size();
	draw_faces_.resize(3*faces_size);
	for (int fi = 0; fi < faces_size; ++fi)
	{
		const Face& f = get_face(fi);
		HalfedgeHandle heh = f.halfedge_handle_;
		Halfedge& he = get_halfedge(heh);
		HalfedgeHandle nheh = he.next_halfedge_handle_;
		HalfedgeHandle pheh = he.prev_halfedge_handle_;
		Halfedge& nhe = get_halfedge(nheh);
		Halfedge& phe = get_halfedge(pheh);
		const IMesh::Point& p0 = get_point(phe.vertex_handle_);
		const IMesh::Point& p1 = get_point( he.vertex_handle_);
		const IMesh::Point& p2 = get_point(nhe.vertex_handle_);
		Normal p1p0 = p0 - p1;
		Normal p1p2 = p2 - p1;
		Normal n = p1p2 % p1p0;
		float norm = n.length();
		Normal unit_normal =
			( (norm != float(0)) ? ((n *= (float(1)/norm)),n) : Normal(0,0,0) );

		vertex_normal_[phe.vertex_handle_] += unit_normal;
		vertex_faces[phe.vertex_handle_] += 1;
		vertex_normal_[ he.vertex_handle_] += unit_normal;
		vertex_faces[ he.vertex_handle_] += 1;
		vertex_normal_[nhe.vertex_handle_] += unit_normal;
		vertex_faces[nhe.vertex_handle_] += 1;

		// Store the faces of this mesh along the way
		draw_faces_[3*fi  ] = phe.vertex_handle_;
		draw_faces_[3*fi+1] =  he.vertex_handle_;
		draw_faces_[3*fi+2] = nhe.vertex_handle_;
	}
	for (int vi = 0; vi < vertices_size; ++vi)
	{
		if (vertex_faces[vi] != 0) vertex_normal_[vi] /= vertex_faces[vi];
	}
}

void ProgressiveMesh::adjust_outgoing_halfedge(VertexHandle vertex_handle)
{
	HalfedgeHandle shh = get_vertex(vertex_handle).halfedge_handle_;
	HalfedgeHandle hh = shh;
	do
	{
		if (!is_valid_face(get_halfedge(hh).face_handle_))
		{
			get_vertex(vertex_handle).halfedge_handle_ = hh;
			break;
		}
		hh = get_opposite_halfedge(hh).next_halfedge_handle_;
	} while (hh != shh);
}

ProgressiveMesh::HalfedgeHandle ProgressiveMesh::build_loop(
	VertexHandle v1, VertexHandle vl, HalfedgeHandle v1vl)
{
	HalfedgeHandle vlv1 = new_edge(vl, v1);
	FaceHandle f1 = new_face();

	HalfedgeHandle o0 = get_opposite_halfedge_handle(v1vl);
	HalfedgeHandle o1 = get_opposite_halfedge_handle(vlv1);
	FaceHandle     f0 = face_handle(v1vl);

	Halfedge&      he0 = get_halfedge(v1vl);
	HalfedgeHandle ph0 = he0.prev_halfedge_handle_;
	HalfedgeHandle nh0 = he0.next_halfedge_handle_;
	Halfedge&     phe0 = get_halfedge(ph0);
	Halfedge&     nhe0 = get_halfedge(nh0);
	Halfedge&      he1 = get_halfedge(vlv1);
	Halfedge&      oe0 = get_opposite_halfedge(v1vl);
	Halfedge&      oe1 = get_opposite_halfedge(vlv1);

	phe0.next_halfedge_handle_ = o1;
	oe1.prev_halfedge_handle_  = ph0;
	oe1.next_halfedge_handle_  = nh0;
	nhe0.prev_halfedge_handle_ = o1;
	he1.next_halfedge_handle_  = v1vl;
	he0.prev_halfedge_handle_  = vlv1;
	he0.next_halfedge_handle_  = vlv1;
	he1.prev_halfedge_handle_  = v1vl;

	oe1.face_handle_ = f0;
	he0.face_handle_ = f1;
	he1.face_handle_ = f1;

	faces_[f1].halfedge_handle_ = v1vl;
	if (f0 != -1) faces_[f0].halfedge_handle_ = o1;

	if(!is_valid_face(oe1.face_handle_))
		get_vertex(v1).halfedge_handle_ = o1; // adjust_outgoing_halfedge of v1;
	if(!is_valid_face(oe0.face_handle_))
		get_vertex(vl).halfedge_handle_ = o0; // adjust_outgoing_halfedge of vl;

	return vlv1;
}

}
