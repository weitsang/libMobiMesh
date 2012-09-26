//
//  TrulyViewDependentPM.cpp
//  MobiMesh
//
//  Created by Ong Yuh Shin on 5/18/11.
//  Copyright 2011 National University of Singapore. All rights reserved.
//

#include <sys/types.h>
#include <sys/mman.h>
#include <sys/fcntl.h>
#include <sys/stat.h>
#include <fstream>
#include <sstream>
#include <errno.h>

#include "MobiMeshConstants.h"
#include "FileReadErrorException.h"
#include "InvalidFileFormatException.h"
#include "InvalidThresholdException.h"
#include "TrulyViewDependentPM.h"

namespace MobiMesh {

TrulyViewDependentPM::TrulyViewDependentPM() :
    need_refined_(false), is_mesh_loaded_(false), mmapped_(false),
    is_active_(false)
{ }

TrulyViewDependentPM::~TrulyViewDependentPM()
{
    unload_mesh();
}

void TrulyViewDependentPM::load_mesh(std::string const &filename)
{
    if (is_mesh_loaded_)
        unload_mesh();

    free_vert_  .reserve (  (MAX_ECVS_PER_FRAME));
    free_face_  .reserve (2*(MAX_ECVS_PER_FRAME));
    free_edge_  .reserve (3*(MAX_ECVS_PER_FRAME));
    free_vplist_.reserve (2*(MAX_ECVS_PER_FRAME));

    freevi = -1;
    freefi = -1;
    freeei = -1;
    freepi = -1;

    mesh_filename_ = filename;

    readVDPM(filename.c_str());
    activate();

    need_refined_ = true;
    is_mesh_loaded_ = true;
}

void TrulyViewDependentPM::unload_mesh()
{
    is_mesh_loaded_ = false;

    free_vert_.clear();
    free_face_.clear();
    free_edge_.clear();
    free_vplist_.clear();

    deactivate();
}

int TrulyViewDependentPM::add_vertex(const Point &_p)
{
    if (freevi == -1) {
        vertices_.push_back(Vertex());
        points_.push_back(_p); 
        size_t nV = n_vertices();
        vertex_params_.resize(nV);
        vertex_status_.resize(nV);
        vertex_weight_.resize(nV);
        resz_index(nV);
        return nV-1;
    } else {
        size_t i0 = free_vert_[freevi--];
        points_[i0] = _p;
        vertex_status_[i0].set_deleted(false);
        return i0;
    }
}

int TrulyViewDependentPM::new_edge(int _start_vh, int _end_vh)
{
    int _heh = n_halfedges();
    Halfedge he0; he0.vertex_handle_ =   _end_vh; hfedges_.push_back(he0);
    Halfedge he1; he1.vertex_handle_ = _start_vh; hfedges_.push_back(he1);
    edge_status_.resize(n_edges());
    return _heh;
}

int TrulyViewDependentPM::new_face()
{
    faces_.push_back(Face());
    face_status_.resize(n_faces());
    return n_faces()-1;
}

int TrulyViewDependentPM::find_halfedge(int _start_vh, int _end_vh)
{
#ifdef DEBUG
    assert(_start_vh != -1 && _end_vh != -1);
#endif

    int shh = vertices_[_start_vh].halfedge_handle_;
    int  hh = shh;
    do {
        if (hh == -1) break;
        if ( hfedges_[hh].vertex_handle_ == _end_vh ) return hh;
        hh = hfedges_[opp(hh)].next_halfedge_handle_;
    } while (hh != shh);

    return -1;
}

void TrulyViewDependentPM::adjust_outgoing_halfedge(int _vh)
{
    int shh = vertices_[_vh].halfedge_handle_;
    int hh = shh;
    do {
        if (hfedges_[hh].face_handle_ == -1)
        {
            vertices_[_vh].halfedge_handle_ = hh;
            break;
        }
        hh = hfedges_[opp(hh)].next_halfedge_handle_;
    } while (hh != shh);
}

int TrulyViewDependentPM::add_face(unsigned int i0, unsigned int i1, unsigned int i2)
{
    int _vertex_handles[3] = { i0, i1, i2 }; uint _vhs_size = 3;

    int                 vh;
    uint                i, ii, n(_vhs_size), vid;
    std::vector<int>    halfedge_handles(n);
    std::vector<bool>   is_new(n), needs_adjust(n, false);
    int                 inner_next, inner_prev,
    outer_next, outer_prev,
    boundary_next, boundary_prev,
    patch_start, patch_end;

    // cache for set_next_halfedge and vertex' set_halfedge
    typedef NextCacheEntryT<int>           NextCacheEntry;
    typedef std::vector<NextCacheEntry>    NextCache;
    NextCache next_cache; next_cache.reserve(3*n);
    assert (n > 2); // don't allow degenerated faces

    // test for topological errors
    for (i=0, ii=1; i<n; ++i, ++ii, ii%=n)
    {
        if ( !is_boundary_vh(_vertex_handles[i]) )
        {
#ifdef DEBUG
            std::cerr << "PolyMeshT::add_face: complex vertex\n";
#endif
            return -1;
        }
        halfedge_handles[i] = find_halfedge(_vertex_handles[i], _vertex_handles[ii]);
        is_new[i] = (halfedge_handles[i] == -1);
        if (!is_new[i] && !is_boundary_heh(halfedge_handles[i]))
        {
#ifdef DEBUG
            std::cerr << "PolyMeshT::add_face: complex edge\n";
#endif
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
                outer_prev = opposite_halfedge_handle(inner_next);
                outer_next = opposite_halfedge_handle(inner_prev);
                boundary_prev = outer_prev;
                do         boundary_prev = opposite_halfedge_handle(next_halfedge_handle(boundary_prev));
                while (!is_boundary_heh(boundary_prev) || boundary_prev==inner_prev);
                boundary_next = next_halfedge_handle(boundary_prev);
                assert(is_boundary_heh(boundary_prev));
                assert(is_boundary_heh(boundary_next));
                // ok ?
                if (boundary_next == inner_next)
                {
#ifdef DEBUG
                    std::cerr << "PolyMeshT::add_face: patch re-linking failed\n";
#endif
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
        if (is_new[i]) halfedge_handles[i] = new_edge(_vertex_handles[i], _vertex_handles[ii]);
    
    // create the face
    int fh = new_face(); set_halfedge_handle_fh(fh, halfedge_handles[n-1]);
    
    // setup halfedges
    for (i=0, ii=1; i<n; ++i, ++ii, ii%=n)
    {
        vh         =  _vertex_handles[ii];
        inner_prev = halfedge_handles[i];
        inner_next = halfedge_handles[ii];
        
        vid = 0;
        if (is_new[i])  vid |= 1;
        if (is_new[ii]) vid |= 2;
        
        if (vid) {
            outer_prev = opposite_halfedge_handle(inner_next);
            outer_next = opposite_halfedge_handle(inner_prev);
            
            // set outer links
            switch (vid) {
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
                    if (halfedge_handle_vh(vh) == -1) {
                        set_halfedge_handle_vh(vh, outer_next);
                        next_cache.push_back(NextCacheEntry(outer_prev, outer_next)); }
                    else {
                        boundary_next = halfedge_handle_vh(vh);
                        boundary_prev = prev_halfedge_handle(boundary_next);
                        next_cache.push_back(NextCacheEntry(boundary_prev, outer_next));
                        next_cache.push_back(NextCacheEntry(outer_prev, boundary_next));}
                    break;
            }
            // set inner link
            next_cache.push_back(NextCacheEntry(inner_prev, inner_next)); }
        else needs_adjust[ii] = (halfedge_handle_vh(vh) == inner_next);
        
        // set face handle
        set_face_handle(halfedge_handles[i], fh);
    }
    // process next halfedge cache
    NextCache::const_iterator ncIt(next_cache.begin()), ncEnd(next_cache.end());
    for (; ncIt != ncEnd; ++ncIt) set_next_halfedge_handle(ncIt->first, ncIt->second);
    // adjust vertices' halfedge handle
    for (i=0; i<n; ++i) if (needs_adjust[i]) adjust_outgoing_halfedge(_vertex_handles[i]);
    
    return fh;
}

void TrulyViewDependentPM::garbage_collection_fast(bool _v, bool _e, bool _f)
{
    int i0, i1, nV(n_vertices()), nE(n_edges()), nF(n_faces());

    // remove deleted vertices
    if (_v && n_vertices() > 0)
    {
        // Note! We need to sort the free vertex pointers before traverse them; (0.5.9)
        sort(free_vert_.begin(), free_vert_.begin()+(freevi+1), &TrulyViewDependentPM::idxcomp);
        i0=0;  i1=nV-1;
        while (1) 
        {
            // Note! We find 1st deleted and last un-deleted, and then swap them; (0.1.0)
            // Note! We mark the free vertices from edge-collapses to speed up garbage_collection; (0.5.9)
            if (freevi == -1)
                break;
            i0 = free_vert_[freevi--];
            while ( vertex_status_[i1].deleted() && i0 < i1 )
                --i1;
            if (i0 >= i1)
                break;

            std::swap(     vertices_[i0],      vertices_[i1]);
            std::swap(       points_[i0],        points_[i1]);
            std::swap(vertex_params_[i0], vertex_params_[i1]);
            std::swap(vertex_status_[i0], vertex_status_[i1]);
            std::swap(vertex_weight_[i0], vertex_weight_[i1]);
            swap_index(              i0 ,                i1 );
            
            vertex_params_[i0]->vh = i0; // Note! We update VertexHandle links since 'garbage_collection()' will shuffle them;

            // update handles of halfedges
            int       sheh = vertices_[i0].halfedge_handle_;
            Halfedge& shehe = hfedges_[sheh];
            int       heh = sheh;
            int       df = -1;
            int       pdf = shehe.face_handle_;
            const OpenMesh::Vec3f& i0n = vertex_params_[i0]->normal;
            do {
                Halfedge&  hehe = hfedges_[heh];
                Halfedge& ohehe = hfedges_[opp(heh)];
                ohehe.vertex_handle_ = i0;
                
                // Note! We embed the updates of GPU-friendly 'draw_faces_' into garbage_collection; (0.5.8)
                // Note! We fix the segmentation fault while updating draw_faces_ on neptune_nonmanifold.vpm; (0.5.9)
                df = ohehe.face_handle_;
                int vd = hehe.vertex_handle_;
                const OpenMesh::Vec3f& vdn = vertex_params_[vd]->normal;
                if ( df >= 0)
                {
                    draw_faces_[ df*3  ] = i0; vpoints_[ df*3  ] = points_[i0]; vnormal_[ df*3  ] = i0n;
                    draw_faces_[ df*3+2] = vd; vpoints_[ df*3+2] = points_[vd]; vnormal_[ df*3+2] = vdn;
                } // ohehe.vertex_handle_ & hehe.vertex_handle_;
                if (pdf >= 0)
                {
                    draw_faces_[pdf*3+1] = vd; vpoints_[pdf*3+1] = points_[vd]; vnormal_[pdf*3+1] = vdn;
                } //  hehe.vertex_handle_;                                               
                
                pdf = df;
                heh = ohehe.next_halfedge_handle_;
            } while (heh != sheh);

            if ( df >= 0)
            {
                int fi =  df*3;
                const OpenMesh::Vec3f& fp0 = vpoints_[fi], fp1 = vpoints_[fi+1], fp2 = vpoints_[fi+2];
                OpenMesh::Vec3f fn = (fp2-fp1) % (fp0-fp1);
                vnormal_[fi] = fn;  vnormal_[fi+1] = fn;  vnormal_[fi+2] = fn;
            }
        };
        nV = vertex_status_[i1].deleted() ? i1 : i1+1;
        vertices_.resize(nV);
        points_.resize(nV);
        vertex_params_.resize(nV);
        vertex_status_.resize(nV);
        vertex_weight_.resize(nV);
        resz_index(nV);
    }
    // remove deleted faces
    if (_f && n_faces() > 0)
    {
        // Note! We need to sort the free face pointers before traverse them; (0.5.9)
        sort(free_face_.begin(), free_face_.begin()+(freefi+1), &TrulyViewDependentPM::idxcomp);
        i0=0;  i1=nF-1;
        while (1)
        {
            // Note! We find 1st deleted and last un-deleted, and then swap them; (0.1.0)
            // Note! We mark the free faces from edge-collapses to speed up garbage_collection; (0.5.9)
            if (freefi == -1) break; i0 = free_face_[freefi--];
            while ( face_status_[i1].deleted() && i0 < i1 )  --i1;
            if (i0 >= i1) break;
            
            // update handles of halfedges
            int   fh = faces_[i1].halfedge_handle_; Halfedge&  fhe = hfedges_[ fh];  fhe.face_handle_ = i0;
            int  pfh =   fhe.prev_halfedge_handle_; Halfedge& pfhe = hfedges_[pfh]; pfhe.face_handle_ = i0;
            int  nfh =   fhe.next_halfedge_handle_; Halfedge& nfhe = hfedges_[nfh]; nfhe.face_handle_ = i0;
            
            std::swap(      faces_[i0],       faces_[i1]);
            std::swap(face_status_[i0], face_status_[i1]);
            
            // Note! We swap the elements of GPU-friendly 'draw_faces_' in garbage_collection; (0.5.8) 
            
            int fv0 = pfhe.vertex_handle_; const OpenMesh::Vec3f& fp0 = points_[fv0];
            int fv1 =  fhe.vertex_handle_; const OpenMesh::Vec3f& fp1 = points_[fv1];
            int fv2 = nfhe.vertex_handle_; const OpenMesh::Vec3f& fp2 = points_[fv2];
            
            int fi = i0*3;
            OpenMesh::Vec3f fn = (fp2-fp1) % (fp0-fp1);
            draw_faces_[fi  ] = fv0; vpoints_[fi  ] = fp0; vnormal_[fi  ] = fn;
            draw_faces_[fi+1] = fv1; vpoints_[fi+1] = fp1; vnormal_[fi+1] = fn;
            draw_faces_[fi+2] = fv2; vpoints_[fi+2] = fp2; vnormal_[fi+2] = fn;
        };
        nF = face_status_[i1].deleted() ? i1 : i1+1;
        faces_.resize(nF);
        face_status_.resize(nF);
        
        draw_faces_.resize(nF*3);
    }
    // remove deleted edges
    if (_e && n_edges() > 0)
    {
        // Note! We need to sort the free edge pointers before traverse them; (0.5.9)
        sort(free_edge_.begin(), free_edge_.begin()+(freeei+1), &TrulyViewDependentPM::idxcomp);
        i0=0; i1=nE-1;
        while (1)
        {
            // Note! We find 1st deleted and last un-deleted, and then swap them; (0.1.0)
            // Note! We mark the free edges from edge-collapses to speed up garbage_collection; (0.5.9)
            if (freeei == -1) break; i0 = free_edge_[freeei--];         
            while ( edge_status_[i1].deleted() && i0 < i1 )  --i1;
            if (i0 >= i1) break;
            
            int     h0 = 2*i0  ; Halfedge&  he0 = hfedges_[i1*2  ];
            int     h1 = 2*i0+1; Halfedge&  he1 = hfedges_[i1*2+1];
            // update handles of vertices
            Vertex& v0 = vertices_[he0.vertex_handle_];
            Vertex& v1 = vertices_[he1.vertex_handle_];
            if (v0.halfedge_handle_ == 2*i1+1) v0.halfedge_handle_ = h1;
            if (v1.halfedge_handle_ == 2*i1  ) v1.halfedge_handle_ = h0;
            // update handles of halfedges (next/prev)
            int    ph0 = he0.prev_halfedge_handle_;
            int    nh0 = he0.next_halfedge_handle_;
            hfedges_[ph0].next_halfedge_handle_ = h0;
            hfedges_[nh0].prev_halfedge_handle_ = h0;
            int    ph1 = he1.prev_halfedge_handle_;
            int    nh1 = he1.next_halfedge_handle_;
            hfedges_[ph1].next_halfedge_handle_ = h1;
            hfedges_[nh1].prev_halfedge_handle_ = h1;
            // update handles of faces
            Face&   f0 = faces_[he0.face_handle_];
            Face&   f1 = faces_[he1.face_handle_];
            if (f0.halfedge_handle_ == 2*i1  ) f0.halfedge_handle_ = h0;
            if (f1.halfedge_handle_ == 2*i1+1) f1.halfedge_handle_ = h1;
            
            std::swap(hfedges_[2*i0  ],  hfedges_[2*i1  ]);
            std::swap(hfedges_[2*i0+1],  hfedges_[2*i1+1]);
            std::swap(edge_status_[i0], edge_status_[i1]);
        };
        nE = edge_status_[i1].deleted() ? i1 : i1+1;
        hfedges_.resize(2*nE);
        edge_status_.resize(nE);
    }

    // Note! We reset all the free memory pointers here just for sure; (0.5.9)
    freevi = -1;
    freeei = -1;
    freefi = -1;
    
    // Note! We need to do garbage_collection on free_vplist_; (0.6.2)
    while (freepi != -1 ) { vplist.erase(free_vplist_[freepi--]); }
}

void TrulyViewDependentPM::vsplit_fast (VPIter vpit, const OpenMesh::Vec3f& p,
                                        const OpenMesh::Vec3f& l_normal,
                                        const OpenMesh::Vec3f& r_normal,
                                        unsigned int l_pos, unsigned int r_pos,
                                        int vl, int vr, int v1vl, int vrv1)
{
    static const VPIter VPEnd = vplist.end();

    int v0 = add_vertex(p), v1 = vpit->vh; ZParam lvp, rvp;

    lvp.vh      = v0;       rvp.vh      = v1;
    lvp.normal  = l_normal; rvp.normal  = r_normal;
    lvp.pos     = l_pos;    rvp.pos     = r_pos;
    lvp.parent  = vpit;     rvp.parent  = vpit;
    lvp.lchild  = VPEnd;    rvp.lchild  = VPEnd;
    lvp.rchild  = VPEnd;    rvp.rchild  = VPEnd;
    
    // Note! We re-allocate the free vparams to vertex-splits; (0.5.9)
    if (freepi == -1) { vpit->lchild = vplist.insert(vplist.end(), lvp);                 }
    else              { vpit->lchild = free_vplist_[freepi--]; (*(vpit->lchild)) = lvp;  }
    if (freepi == -1) { vpit->rchild = vplist.insert(vplist.end(), rvp);                 }
    else              { vpit->rchild = free_vplist_[freepi--]; (*(vpit->rchild)) = rvp;  }
    
    // vertex_split_fast(v0, v1, vl, vr, v1vl, vrv1);
    
    int vlv0 = -1, v0vr = -1, v0v1 = -1, v1v0 = -1;
    int  f01 = -1,  f10 = -1; 
    
    // Note! The routine finds the halfedges (v1, vl), (vr, v1) in one neighborhood traverse; (0.5.4)
    // Note! We combine the routine with 'get_active_cuts()' to save one neighborhood traverse; (0.5.5)
    
    // build loop from halfedge v1->vl, vlv0 = insert_loop(v1vl);
    if (vl != -1)
    {
        // assert(v1vl != -1);
        
        // Note! We re-allocate the free edge to vertex-splits before request a new one; (0.5.9)
        vlv0 = -1; if (freeei == -1) { vlv0 = n_halfedges(); hfedges_.resize(vlv0+2);           }
        else { vlv0 = free_edge_[freeei--]; edge_status_[vlv0].set_deleted(false); vlv0 <<= 1;  } 
        
        int         v0vl = vlv0+1;
        Halfedge&   hel0 = hfedges_[vlv0]; hel0.vertex_handle_ = v1;
        Halfedge&   he0l = hfedges_[v0vl]; he0l.vertex_handle_ = vl;           
        
        // Note! We re-allocate the free face to vertex-splits before request a new one; (0.5.9)
        if (freefi == -1) { faces_.push_back(Face()); draw_faces_.resize(n_faces()*3); f01 = n_faces()-1; } 
        else { f01 = free_face_[freefi--]; face_status_[f01].set_deleted(false); }
        
        int         vlv1 = (v1vl&1) ? v1vl-1 : v1vl+1;
        Halfedge&   hel1 = hfedges_[vlv1];
        Halfedge&   he1l = hfedges_[v1vl];
        int          f1l = he1l.face_handle_;
        int        phh1l = he1l.prev_halfedge_handle_;
        int        nhh1l = he1l.next_halfedge_handle_;
        Halfedge&  phe1l = hfedges_[phh1l];
        Halfedge&  nhe1l = hfedges_[nhh1l];
        
        phe1l.next_halfedge_handle_ = v0vl; 
        he0l.prev_halfedge_handle_ = phh1l;
        he0l.next_halfedge_handle_ = nhh1l; 
        nhe1l.prev_halfedge_handle_ = v0vl;
        hel0.next_halfedge_handle_ = v1vl;
        he1l.prev_halfedge_handle_ = vlv0;
        he1l.next_halfedge_handle_ = vlv0;
        hel0.prev_halfedge_handle_ = v1vl;
        
        he0l.face_handle_ = f1l;
        he1l.face_handle_ = f01;
        hel0.face_handle_ = f01;
        
        faces_[f01].halfedge_handle_ = v1vl;
        if (f1l != -1) faces_[f1l].halfedge_handle_ = v0vl;
        
        if(he0l.face_handle_ == -1) vertices_[v1].halfedge_handle_ = v0vl; // adjust_outgoing_halfedge of v1;
        if(hel1.face_handle_ == -1) vertices_[vl].halfedge_handle_ = vlv1; // adjust_outgoing_halfedge of vl;
    }
    
    // build loop from halfedge vr->v1, insert_loop(vrv1);
    if (vr != -1)
    {
        // assert(vrv1 != -1);
        
        // Note! We re-allocate the free edge to vertex-splits before request a new one; (0.5.9)
        v0vr = -1; if (freeei == -1) { v0vr = n_halfedges(); hfedges_.resize(v0vr+2);           }
        else { v0vr = free_edge_[freeei--]; edge_status_[v0vr].set_deleted(false); v0vr <<= 1;  }
        
        int         vrv0 = v0vr+1;
        Halfedge&   he0r = hfedges_[v0vr]; he0r.vertex_handle_ = vr;
        Halfedge&   her0 = hfedges_[vrv0]; her0.vertex_handle_ = v1;
        
        // Note! We re-allocate the free face to vertex-splits before request a new one; (0.5.9)
        if (freefi == -1) { faces_.push_back(Face()); draw_faces_.resize(n_faces()*3); f10 = n_faces()-1; } 
        else { f10 = free_face_[freefi--]; face_status_[f10].set_deleted(false); }
        
        int         v1vr = (vrv1&1) ? vrv1-1 : vrv1+1;
        Halfedge&   he1r = hfedges_[v1vr];
        Halfedge&   her1 = hfedges_[vrv1];
        int          fr1 = her1.face_handle_;
        int        phhr1 = her1.prev_halfedge_handle_;
        int        nhhr1 = her1.next_halfedge_handle_;
        Halfedge&  pher1 = hfedges_[phhr1]; 
        Halfedge&  nher1 = hfedges_[nhhr1];
        
        pher1.next_halfedge_handle_ = vrv0;  
        her0.prev_halfedge_handle_ = phhr1;
        her0.next_halfedge_handle_ = nhhr1;  
        nher1.prev_halfedge_handle_ = vrv0;
        he0r.next_halfedge_handle_ = vrv1;
        her1.prev_halfedge_handle_ = v0vr;
        her1.next_halfedge_handle_ = v0vr;
        he0r.prev_halfedge_handle_ = vrv1;
        
        her0.face_handle_ = fr1;
        her1.face_handle_ = f10;
        he0r.face_handle_ = f10;
        
        faces_[f10].halfedge_handle_ = vrv1;
        if (fr1 != -1) faces_[fr1].halfedge_handle_ = vrv0;
        
        if(her0.face_handle_ == -1) vertices_[vr].halfedge_handle_ = vrv0; // adjust_outgoing_halfedge of vr;
        if(he1r.face_handle_ == -1) vertices_[v1].halfedge_handle_ = v1vr; // adjust_outgoing_halfedge of v1;
    }
    
    // handle boundary cases
    if (vl == -1) vlv0 = prev_halfedge_handle(halfedge_handle_vh(v1));
    if (vr == -1) vrv1 = prev_halfedge_handle(halfedge_handle_vh(v1));
    // assert(vlv0 != -1 && vrv1 != -1 && v1 == to_vertex_handle(vrv1));
    
    // split vertex v1 into edge v0v1, v0v1 = insert_edge(v0, vlv0, vrv1);
    
    // Note! We re-allocate the free edge to vertex-splits before request a new one; (0.5.9)
    v0v1 = -1; if (freeei == -1) { v0v1 = n_halfedges(); hfedges_.resize(v0v1+2);           }
    else { v0v1 = free_edge_[freeei--]; edge_status_[v0v1].set_deleted(false); v0v1 <<= 1;  }
    
    v1v0 = v0v1+1;
    Halfedge&       he01 = hfedges_[v0v1]; he01.vertex_handle_ = v1;
    Halfedge&       he10 = hfedges_[v1v0]; he10.vertex_handle_ = v0;
    
    Halfedge&       hel0 = hfedges_[vlv0]; 
    Halfedge&       her1 = hfedges_[vrv1];
    int             nhl1 = hel0.next_halfedge_handle_;
    int             nhr1 = her1.next_halfedge_handle_;
    Halfedge&      nhel1 = hfedges_[nhl1];
    Halfedge&      nher1 = hfedges_[nhr1]; 
    
    he01.next_halfedge_handle_ = nhl1;
    nhel1.prev_halfedge_handle_ = v0v1;
    hel0.next_halfedge_handle_ = v0v1;
    he01.prev_halfedge_handle_ = vlv0;
    he10.next_halfedge_handle_ = nhr1;
    nher1.prev_halfedge_handle_ = v1v0;
    her1.next_halfedge_handle_ = v1v0;
    he10.prev_halfedge_handle_ = vrv1;

    he01.face_handle_ = f01; if (f01 != -1) faces_[f01].halfedge_handle_ = v0v1;
    he10.face_handle_ = f10; if (f10 != -1) faces_[f10].halfedge_handle_ = v1v0;

    // we do not always need to adjust the outgoing (boundary) halfedge of v1.
    Vertex& vv0 = vertices_[v0];
    Vertex& vv1 = vertices_[v1];
    const OpenMesh::Vec3f& v0n = l_normal;

    int   hehv1 = vv1.halfedge_handle_; vv0.halfedge_handle_ = v0v1;
    int     heh = v0v1; int df = -1, pdf = he01.face_handle_; 
    do
    {
        Halfedge&   hehe = hfedges_[   (heh)]; 
        Halfedge&  ohehe = hfedges_[opp(heh)];
        if (heh == hehv1)            {vv1.halfedge_handle_ = v1v0;} // adjust outgoing halfedge of v1
        if (hehe.face_handle_ == -1) {vv0.halfedge_handle_ = heh ;} // adjust outgoing halfedge of v0
        ohehe.vertex_handle_ = v0;
        
        // Note! We embed the updates of GPU-friendly 'draw_faces_' into vertex-splits; (0.5.8)
        // Note! We fix the segmentation fault while updating draw_faces_ on neptune_nonmanifold.vpm; (0.5.9)
        df = ohehe.face_handle_;
        int vd = hehe.vertex_handle_;
        const OpenMesh::Vec3f& vdn = vertex_params_[vd]->normal;
        if (df >= 0)
        {
            draw_faces_[ df*3  ] = v0; vpoints_[ df*3  ] = points_[v0]; vnormal_[ df*3  ] = v0n;
            draw_faces_[ df*3+2] = vd; vpoints_[ df*3+2] = points_[vd]; vnormal_[ df*3+2] = vdn;
        } // ohehe.vertex_handle_ & hehe.vertex_handle_;
        if (pdf >= 0)
        {
            draw_faces_[pdf*3+1] = vd; vpoints_[pdf*3+1] = points_[vd]; vnormal_[pdf*3+1] = vdn;
        } //  hehe.vertex_handle_; 
        pdf = df;
        heh = ohehe.next_halfedge_handle_;
    } while (heh != v0v1);
    
    if ( df >= 0)
    {
        int fi =  df*3;
        const OpenMesh::Vec3f& fp0 = vpoints_[fi], fp1 = vpoints_[fi+1], fp2 = vpoints_[fi+2]; 
        OpenMesh::Vec3f fn = (fp2-fp1) % (fp0-fp1);
        vnormal_[fi]   = fn;
        vnormal_[fi+1] = fn;
        vnormal_[fi+2] = fn;
    }

    // make sure the outgoing halfedge adjustment conforms to the assumption
    // (only one boundary per vertex)    
    vertex_params_[v0] = vpit->lchild;
    vertex_params_[v1] = vpit->rchild;

    // split the original vertex weight between the vertices
    vertex_weight_[v0] = vertex_weight_[v1]/2;
    vertex_weight_[v1] = vertex_weight_[v1]/2;

    splt_index(v0,v1);
}

bool TrulyViewDependentPM::ecol_fast (VPIter vpit)
{
    static const VPIter VPEnd = vplist.end();

    int   v0 = vpit->lchild->vh;
    int   v1 = vpit->rchild->vh;
    int v0v1 = -1, vl = -1, vr = -1;

    if (vertex_status_[v0].deleted() || vertex_status_[v1].deleted())
        return false;

    Vertex& vv0 = vertices_[v0]; 
    int v0hh = vv0.halfedge_handle_;

#ifdef DEBUG
    assert(v0hh != -1);
#endif

    // Find halfedge (v0, v1) here and remark v0's neighbors before
    // intersection tests; (0.5.6)
    ni = 0; int heh = v0hh, vh = -1;
    do
    {
        vh = hfedges_[heh].vertex_handle_;
        vertex_status_[vh].set_tagged(false);
        neighbors[ni++] = vh;
        if (vh == v1)
            v0v1 = heh;
        heh = hfedges_[opp(heh)].next_halfedge_handle_;
    } while (heh != v0hh);

#ifdef DEBUG
    assert (ni < MAX_NEIGHBORS);
#endif

    int         v1v0 = opp(v0v1);
    Halfedge&   he01 = hfedges_[v0v1];
    Halfedge&   he10 = hfedges_[v1v0];

    int         v1vl = he01.next_halfedge_handle_;  
    int         vlv0 = he01.prev_halfedge_handle_;  
    int         vlv1 = opp(v1vl);
    int         v0vl = opp(vlv0);
    Halfedge&   hel1 = hfedges_[vlv1]; 
    Halfedge&   he0l = hfedges_[v0vl]; 

    // the edges v1-vl and vl-v0 must not be both boundary edges
    if (he01.face_handle_ != -1)
    {
        if (he0l.face_handle_ == -1 && hel1.face_handle_ == -1)
            return false;
        vl = he0l.vertex_handle_;
    }

    int         v0vr = he10.next_halfedge_handle_;  
    int         vrv1 = he10.prev_halfedge_handle_;  
    int         vrv0 = opp(v0vr);
    int         v1vr = opp(vrv1);
    Halfedge&   her0 = hfedges_[vrv0];
    Halfedge&   he1r = hfedges_[v1vr];

    // the edges v0-vr and vr-v1 must not be both boundary edges
    if (he10.face_handle_ != -1)
    {
        if (he1r.face_handle_ == -1 && her0.face_handle_ == -1)
            return false;
        vr = he1r.vertex_handle_;
    }

    // if vl and vr are equal or both invalid -> fail
    if (vl == vr)
        return false;

    Vertex&      vv1 = vertices_[v1]; 
    const int&  v1hh = vv1.halfedge_handle_;
#ifdef DEBUG
    assert(v1hh!=-1);
#endif
    Halfedge&   v0he = hfedges_[v0hh];
    Halfedge&   v1he = hfedges_[v1hh];
    
    // edge between two boundary vertices should be a boundary edge
    if (   v0he.face_handle_ == -1 // is_boundary_vh(v0)
        && v1he.face_handle_ == -1 // is_boundary_vh(v1)
        && he01.face_handle_ != -1 
        && he10.face_handle_ != -1)
        return false;

    heh = v1hh;
    do
    {
        vh = hfedges_[heh].vertex_handle_;
        vertex_status_[vh].set_tagged(true);
        heh = hfedges_[opp(heh)].next_halfedge_handle_;
    } while (heh != v1hh);
    
    for (size_t i = 0; i < ni; i++)
    {
        vh = neighbors[i];
        if (vertex_status_[vh].tagged() && vh != vl && vh != vr)
            return false;
    }
    
    Halfedge&   he1l = hfedges_[v1vl]; 
    Halfedge&   hel0 = hfedges_[vlv0];
    Halfedge&   he0r = hfedges_[v0vr];
    Halfedge&   her1 = hfedges_[vrv1];
    
    int   f01 = he01.face_handle_; 
    int   f10 = he10.face_handle_; 

    heh = v0v1; int df = -1, pdf = he01.face_handle_;
    const OpenMesh::Vec3f& v1n = vpit->normal;
    do {
        Halfedge&  hehe = hfedges_[   (heh)];
        Halfedge& ohehe = hfedges_[opp(heh)];
        ohehe.vertex_handle_ = v1;

        // Note! We embed the updates of GPU-friendly 'draw_faces_' into edge-collapses; (0.5.8)
        // Note! We fix the segmentation fault while updating draw_faces_ on neptune_nonmanifold.vpm; (0.5.9)
        df = ohehe.face_handle_;
        int vd = hehe.vertex_handle_;
        const OpenMesh::Vec3f& vdn = vertex_params_[vd]->normal;
        if ( df >= 0)
        {
            draw_faces_[ df*3  ] = v1; vpoints_[ df*3  ] = points_[v1]; vnormal_[ df*3  ] = v1n;
            draw_faces_[ df*3+2] = vd; vpoints_[ df*3+2] = points_[vd]; vnormal_[ df*3+2] = vdn;
        } // ohehe.vertex_handle_ & hehe.vertex_handle_;
        if (pdf >= 0)
        {
            draw_faces_[pdf*3+1] = vd; vpoints_[pdf*3+1] = points_[vd]; vnormal_[pdf*3+1] = vdn;
        } //  hehe.vertex_handle_; 
        
        pdf = df;
        heh = ohehe.next_halfedge_handle_;
    } while (heh != v0v1);

    // adjust outgoing halfedge of v1
    if (v0hh != v0v1 && v0he.face_handle_ == -1)
        vv1.halfedge_handle_ = v0hh;
    if (v1hh == v1v0)
        vv1.halfedge_handle_ = v1vl;
    
    hel0.next_halfedge_handle_ = v1vl;
    he1l.prev_halfedge_handle_ = vlv0;
    her1.next_halfedge_handle_ = v0vr;
    he0r.prev_halfedge_handle_ = vrv1;
    
    if (f01 != -1)
        faces_[f01].halfedge_handle_ = v1vl;
    if (f10 != -1)
        faces_[f10].halfedge_handle_ = v0vr;
    
    // Note! We mark the free edge and vertex from edge-collapses; (0.5.9)
    vv0.halfedge_handle_ = -1;
    edge_status_[v0v1>>1].set_deleted(true);
    free_edge_[++freeei] = (v0v1>>1);
    vertex_status_[v0].set_deleted(true);
    free_vert_[++freevi] = v0;
    vertex_params_[v0] = VPEnd;

    // remove loops, collapse_loop(next_halfedge_handle(hn)) 
    if (he1l.next_halfedge_handle_ == vlv0
        && hel0.prev_halfedge_handle_ == v1vl) 
    {
        int pv0vl = he0l.prev_halfedge_handle_;
        Halfedge& phe0l = hfedges_[pv0vl];
        int nv0vl = he0l.next_halfedge_handle_;
        Halfedge& nhe0l = hfedges_[nv0vl];

        int  f0l = he0l.face_handle_;

        he1l.next_halfedge_handle_ = nv0vl;
        nhe0l.prev_halfedge_handle_ =  v1vl;
        he1l.prev_halfedge_handle_ = pv0vl;
        phe0l.next_halfedge_handle_ =  v1vl;

        he1l.face_handle_ = f0l;

        Vertex& vvl = vertices_[vl];
        if (vv1.halfedge_handle_ == v0vl)
            vv1.halfedge_handle_ =  v1vl; // adjust outgoing halfedge of v1
        if (vvl.halfedge_handle_ == vlv0)
            vvl.halfedge_handle_ = nv0vl; // adjust outgoing halfedge of vl

        if (f0l != -1 && faces_[f0l].halfedge_handle_ == v0vl)
            faces_[f0l].halfedge_handle_ = v1vl;

        // Note! We mark the free face and edge from edge-collapses; (0.5.9)
        if (f01 != -1)
        {
            faces_[f01].halfedge_handle_ = -1;
            face_status_[f01].set_deleted(true);
            free_face_[++freefi] = f01;
        }
        edge_status_[vlv0>>1].set_deleted(true);
        free_edge_[++freeei] = (vlv0>>1);
    }

    // remove loops, collapse_loop(on)
    if (he0r.next_halfedge_handle_ == vrv1
        && her1.prev_halfedge_handle_ == v0vr) 
    {
        int  pvrv0 = her0.prev_halfedge_handle_;
        Halfedge& pher0 = hfedges_[pvrv0];
        int  nvrv0 = her0.next_halfedge_handle_;
        Halfedge& nher0 = hfedges_[nvrv0];
        
        int  fr0 = her0.face_handle_;
        
        her1.next_halfedge_handle_ = nvrv0;
        nher0.prev_halfedge_handle_ =  vrv1;
        her1.prev_halfedge_handle_ = pvrv0;
        pher0.next_halfedge_handle_ =  vrv1;
        
        her1.face_handle_ = fr0;

        Vertex& vvr = vertices_[vr];
        if (vvr.halfedge_handle_ == vrv0)
            vvr.halfedge_handle_ =  vrv1; // adjust outgoing halfedge of vr
        if (vv1.halfedge_handle_ == v0vr)
            vv1.halfedge_handle_ = nvrv0; // adjust outgoing halfedge of v1
        
        if (fr0 != -1 && faces_[fr0].halfedge_handle_ == vrv0)
            faces_[fr0].halfedge_handle_ = vrv1;

        // Note! We mark the free face and edge from edge-collapses; (0.5.9)
        if (f10 != -1)
        {
            faces_[f10].halfedge_handle_ = -1;
            face_status_[f10].set_deleted(true);
            free_face_[++freefi] = f10;
        }
        edge_status_[v0vr>>1].set_deleted(true);
        free_edge_[++freeei] = (v0vr>>1);
    }

    vertex_params_[v1] = vpit;
    vpit->vh = v1;
    vertex_weight_[v1] = vertex_weight_[v0] + vertex_weight_[v1];

    // Note! We mark the free vparams from edge-collapses; (0.5.9)
    free_vplist_[++freepi] = vpit->lchild;
    vpit->lchild = VPEnd;
    free_vplist_[++freepi] = vpit->rchild;
    vpit->rchild = VPEnd;
    
    ecol_index(v0,v1);

    return true;
}

bool TrulyViewDependentPM::adaptive_refinement(int max_refinements,
                                               unsigned int edge_collapse_threshold,
                                               unsigned int vertex_split_threshold)
{
    // will be set to true if we hit the max number of refinements for either
    // the vertex splits or edge collapses (i.e. means there is more to refine
    // despite this function having returned)
    bool need_further_refined_ = false;

    static const VPIter VPEnd = vplist.end();
    
    // edge-collapse (not root)
    size_t eccount = 0;
    for (size_t vit = 0; vit != n_vertices(); )
    {
        if(eccount >= max_refinements)
        {
            need_further_refined_ = true;
            break;
        }

        VPIter vpit = vertex_params_[vit];

        if ( vpit != VPEnd ) 
        {
            VPIter pvpit = vpit->parent;

            // Only query the parent from rchild (not lchild) to save some
            // visual queries
            if (pvpit != VPEnd 
                && pvpit->rchild == vpit 
                && pvpit->lchild->lchild == VPEnd)
            {
                pvpit->vh = vit;
                // Check that both vertices involved in the collapse have low
                // enough vertex weight.
                int left_vh = pvpit->lchild->vh;
                int right_vh = pvpit->rchild->vh;
                if (vertex_weight_[left_vh] >= edge_collapse_threshold
                    || vertex_weight_[right_vh] >= edge_collapse_threshold)
                {
                    ++vit;
                    continue;
                }

                if (!ecol_fast (pvpit))
                    ++vit;
                else
                    ++eccount;
            }
            else
                ++vit;
        }
        else
            ++vit;
    }

    // vertex-split (not leaf)
    size_t vscount = 0;
    for (size_t vit = 0; vit != n_vertices(); )
    {
        if (vertex_weight_[vit] <= vertex_split_threshold)
        {
            ++vit;
            continue;
        }

        if (vscount >= max_refinements)
        {
            need_further_refined_ = true;
            break;
        }

        VPIter vpit = vertex_params_[vit]; 
        if ( vpit != VPEnd && vpit->pos != 0)
        {
            force_vsplit(vpit);
            ++vscount;
        }
        else
            ++vit;
    }

    edge_status_.resize(n_edges());
    face_status_.resize(n_faces());

    garbage_collection_fast(true, true, true);

    return need_further_refined_;
}

void TrulyViewDependentPM::activate()
{
    if (is_active())
        return;

    // Get the mmap size
    struct stat sbuf;
    if (stat(mesh_filename_.c_str(), &sbuf) == -1)
		throw FileReadErrorException("stat() read error");
    map_size = sbuf.st_size;

    // Open the mesh file for reading
    if ((fd = open(mesh_filename_.c_str(), O_RDONLY)) == -1)
		throw FileReadErrorException("open() failed, unable to open file");

    // Use mmap() if possible so that we can access the file easily
    if ((map_data = (char *) mmap(NULL, map_size, PROT_READ, MAP_SHARED, fd, 0))
        != MAP_FAILED)
	{
        mmapped_ = true;
#ifdef DEBUG
        std::cout << "mmap() succeeded" << std::endl;
#endif
	}
    else
    {
        mmapped_ = false;
        map_data = NULL;
#ifdef DEBUG
        std::cout << "mmap() failed, vertices will be loaded using pread() "
                     "instead" << std::endl;
#endif
    }

    is_active_ = true;
}

void TrulyViewDependentPM::deactivate()
{
    if (!is_active())
        return;

    is_active_ = false;

    if (mmapped_)
    {
        if (munmap(map_data, map_size) == -1)
        {
#ifdef DEBUG
            std::cerr << "Error un-mmapping" << std::endl;
#endif
        }
        mmapped_ = false;
        map_data = NULL;
    }

    close(fd);
}

bool TrulyViewDependentPM::adapt(int max_refinements,
                                 unsigned int edge_collapse_threshold,
                                 unsigned int vertex_split_threshold)
{
    if (2*edge_collapse_threshold >= vertex_split_threshold)
    {
        throw InvalidThresholdException("edge collapse threshold must be "
                                        "strictly less than half of vertex "
                                        "split threshold");
    }

    max_refinements = (max_refinements > MAX_ECVS_PER_FRAME) ?
        MAX_ECVS_PER_FRAME : max_refinements;

    need_refined_ = adaptive_refinement(max_refinements,
                                        edge_collapse_threshold,
                                        vertex_split_threshold);

    return need_refined_;
}

}