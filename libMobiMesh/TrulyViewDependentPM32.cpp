//
//  TrulyViewDependentPM32.cpp
//  MobiMesh
//
//  Created by Ong Yuh Shin on 5/18/11.
//  Copyright 2011 National University of Singapore. All rights reserved.
//

#include <sys/mman.h>
#include <sys/fcntl.h>
#include <fstream>
#include <sstream>
#include <errno.h>

#include "TrulyViewDependentPM32.h"
#include "FileReadErrorException.h"
#include "InvalidFileFormatException.h"

namespace MobiMesh {

TrulyViewDependentPM32::TrulyViewDependentPM32() { }

TrulyViewDependentPM32::TrulyViewDependentPM32(const std::string& filename)
    : TrulyViewDependentPM()
{
    load_mesh(filename);
}

TrulyViewDependentPM32::~TrulyViewDependentPM32() { }

void TrulyViewDependentPM32::get_active_cuts_fast(VPIter vpit,
                                                  unsigned int fund_lcut_tree,
                                                  unsigned int fund_rcut_tree,
                                                  unsigned int fund_lcut_node,
                                                  unsigned int fund_rcut_node,
                                                  int &vl, int &vr,
                                                  int &v1vl, int &vrv1)
{
    // Note! We remove the dependencies on VHierarchyNodeIndex to save some CPU time; (0.5.5)
    vl = -1; vr = -1; v1vl = -1; vrv1 = -1;
    
    int sheh = vertices_[vpit->vh].halfedge_handle_;
    int  heh = sheh;
    do
    {
        const Halfedge&        hehe = hfedges_[heh]; 
        int                   vv_it = hehe.vertex_handle_;

        VertexID           neighbor = vertex_vindex_[vv_it];
        unsigned int  neighbor_tree = neighbor >> tree_id_bitshift_;
        unsigned int  neighbor_node = neighbor &  tree_id_mask_    ; 

        // Note! We finds the halfedges (v1, vl), (vr, v1) here to save one neighborhood traverse; (0.5.5)

        if (vl == -1 && neighbor_tree == fund_lcut_tree)
        {
            VertexID lcid = fund_lcut_node;
            while (lcid > neighbor_node) { lcid >>= 1; } 
            if (neighbor_node == lcid)
            {
                vl = vv_it; v1vl = heh;
            } // v1vl = find_halfedge(v1, vl);
        }
        if (vr == -1 && neighbor_tree == fund_rcut_tree)
        {
            VertexID rcid = fund_rcut_node;
            while (rcid > neighbor_node) { rcid >>= 1; }
            if (neighbor_node == rcid) 
            {
                vr = vv_it;
                vrv1 = opposite_halfedge_handle(heh);
            } // vrv1 = find_halfedge(vr, v1);
        }
        if (vl != -1 && vr != -1)
            break;

        heh = hfedges_[opp(heh)].next_halfedge_handle_;
    } while (heh != sheh);
}

TrulyViewDependentPM32::ZsData TrulyViewDependentPM32::loadVs(unsigned int pos)
{
    if (mmapped_)
        return (ZsData &) * (map_data + (s_base_ + (pos-1)*s_details_));
    else
    {
        ZsData data;
        pread(fd, &data, sizeof(ZsData), (s_base_ + (pos-1)*s_details_));
        return data;
    }
}

void TrulyViewDependentPM32::force_vsplit (VPIter vpit)
{
    const ZsData& zs = loadVs(vpit->pos);

    unsigned int    fund_lcut      = zs.fund_lcut_index;
    unsigned int    fund_rcut      = zs.fund_rcut_index;
    unsigned int    fund_lcut_tree = fund_lcut >> (tree_id_bitshift_);
    unsigned int    fund_rcut_tree = fund_rcut >> (tree_id_bitshift_);
    unsigned int    fund_lcut_node = fund_lcut &  (tree_id_mask_);
    unsigned int    fund_rcut_node = fund_rcut &  (tree_id_mask_);

    int vl = -1, vr = -1, v1vl = -1, vrv1 = -1;
    get_active_cuts_fast (vpit,fund_lcut_tree,fund_rcut_tree,fund_lcut_node,fund_rcut_node,vl,vr,v1vl,vrv1);
    
    while (vl == vr) {
        force_vsplit (vertex_params_[vl]);
        get_active_cuts_fast (vpit,fund_lcut_tree,fund_rcut_tree,fund_lcut_node,fund_rcut_node,vl,vr,v1vl,vrv1);
    }

    int pos = zs.c_pos; unsigned int l_pos = 0, r_pos = 0;
    if      (pos == 0    ) { l_pos = 0;                       r_pos = 0;                        }
    else if (pos  < 0    ) { l_pos = 0;                       r_pos = (unsigned int) -pos >> 1; }
    else if (pos % 2 == 1) { l_pos = (unsigned int) pos >> 1; r_pos = l_pos + 1;                }
    else                   { l_pos = (unsigned int) pos >> 1; r_pos = 0;                        }

    vsplit_fast(vpit,zs.p,zs.l_normal,zs.r_normal,l_pos,r_pos,vl,vr,v1vl,vrv1);
}

inline void TrulyViewDependentPM32::splt_index (int v0, int v1)
{
    VertexID    parent_id = vertex_vindex_[v1];
    VertexID      tree_id = parent_id >> tree_id_bitshift_;
    VertexID      node_id = parent_id  & tree_id_mask_    ;
    VertexID    lchild_id = ( tree_id << tree_id_bitshift_) | ((node_id<<1)  );
    VertexID    rchild_id = ( tree_id << tree_id_bitshift_) | ((node_id<<1)+1);
    
    vertex_vindex_[v0] = lchild_id;
    vertex_vindex_[v1] = rchild_id;
}

inline void TrulyViewDependentPM32::ecol_index (int v0, int v1) 
{
    VertexID    rchild_id = vertex_vindex_[v1];
    VertexID      tree_id = rchild_id >> tree_id_bitshift_;
    VertexID      node_id = rchild_id  & tree_id_mask_    ;
    VertexID    parent_id = ( tree_id << tree_id_bitshift_) | (node_id>>1);
    
    vertex_vindex_[v0] = 0;
    vertex_vindex_[v1] = parent_id;
}

void TrulyViewDependentPM32::readVDPM(const char* _filename)
{
    char                        fileformat[16];
    unsigned int                fvi[3], pos;
    OpenMesh::Vec3f             p, normal;
    
    vpoints_.reserve(MAX_BUFF_FACES*3);
    vpoints_.resize (MAX_BUFF_FACES*3);
    vnormal_.reserve(MAX_BUFF_FACES*3);
    vnormal_.resize (MAX_BUFF_FACES*3);
    
    std::ifstream ifs(_filename, std::ios::binary);
    if (!ifs)
		throw FileReadErrorException("ifstream read error");

    ifs.read(fileformat, 10);
    fileformat[10] = '\0';
    if (std::string(fileformat) != std::string("VDProgMesh"))
    {
		ifs.close();
		throw InvalidFileFormatException("Invalid file format", "VDProgMesh", fileformat);
    }

    ifs.read( (char*)&n_base_vertices_, sizeof(n_base_vertices_) );
    ifs.read( (char*)&n_base_faces_   , sizeof(n_base_faces_) );
    ifs.read( (char*)&n_details_      , sizeof(n_details_) );

    calc_tree_id_bits (n_base_vertices_);
    for (size_t i=0; i<n_base_vertices_; ++i)
    {
        ifs.read( (char*)&p            , sizeof(p) );
        ifs.read( (char*)&normal       , sizeof(normal) );
        ifs.read( (char*)&pos          , sizeof(pos) );
        
        ZParam vp;
        vp.vh            = add_vertex(p);
        vp.normal        = normal;
        vp.pos           = pos;
        vp.parent        = vplist.end(); 
        vp.lchild        = vplist.end(); 
        vp.rchild        = vplist.end();
        
        vertex_params_[vp.vh] = vplist.insert(vplist.end(),vp);
        vertex_vindex_[vp.vh] = ( (VertexID) i << tree_id_bitshift_) | 1; 
        vertex_weight_[vp.vh] = 0;
    }
    draw_faces_.resize(n_base_faces_*3);
    for (size_t i=0; i<n_base_faces_; ++i)
    {
        ifs.read( (char*)&fvi[0] , sizeof(fvi[0]) );
        ifs.read( (char*)&fvi[1] , sizeof(fvi[1]) );
        ifs.read( (char*)&fvi[2] , sizeof(fvi[2]) );
        
        add_face(fvi[0],fvi[1],fvi[2]);
        
        int fv0 = fvi[0]; const OpenMesh::Vec3f& fp0 = points_[fv0]; const OpenMesh::Vec3f& fn0 = vertex_params_[fv0]->normal; 
        int fv1 = fvi[1]; const OpenMesh::Vec3f& fp1 = points_[fv1]; const OpenMesh::Vec3f& fn1 = vertex_params_[fv1]->normal;
        int fv2 = fvi[2]; const OpenMesh::Vec3f& fp2 = points_[fv2]; const OpenMesh::Vec3f& fn2 = vertex_params_[fv2]->normal;
        
        int fi = i*3; 
        draw_faces_[fi  ] = fv0; vpoints_[fi  ] = fp0; vnormal_[fi  ] = fn0;
        draw_faces_[fi+1] = fv1; vpoints_[fi+1] = fp1; vnormal_[fi+1] = fn1;
        draw_faces_[fi+2] = fv2; vpoints_[fi+2] = fp2; vnormal_[fi+2] = fn2;
    }
    ifs.close();
    
    // Note! We use vertex normals from .vpm directly without real-time updates; (0.5.2)
    // update_face_normals();
    
    // Note! We calculate the base offset and data size to access mapped files properly; (0.5.7)
    
    size_t s_params_, s_base_vertices_, s_base_faces_;
    
    s_params_        = sizeof("VDProgMesh")-1+sizeof(unsigned int)+sizeof(unsigned int)+sizeof(unsigned int); 
    s_base_vertices_ = sizeof(OpenMesh::Vec3f)+sizeof(OpenMesh::Vec3f)+sizeof(unsigned int);
    s_base_faces_    = sizeof(unsigned int)+sizeof(unsigned int)+sizeof(unsigned int);
    
    s_base_          = s_params_+n_base_vertices_*s_base_vertices_+n_base_faces_*s_base_faces_;
    s_details_       = sizeof(OpenMesh::Vec3f)+sizeof(OpenMesh::Vec3f)+sizeof(OpenMesh::Vec3f)+sizeof(unsigned int)+sizeof(unsigned int)+sizeof(int);

#ifdef DEBUG
    assert (s_details_ == sizeof(ZsData)); 
#endif
}

}