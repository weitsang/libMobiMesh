//
//  TrulyViewDependentPM64.h
//  MobiMesh
//
//  Created by Ong Yuh Shin on 5/18/11.
//  Copyright 2011 National University of Singapore. All rights reserved.
//

#ifndef __TRULY_VIEW_DEPENDENT_PM64_H_
#define __TRULY_VIEW_DEPENDENT_PM64_H_

#include <sys/stat.h>

#include "TrulyViewDependentPM.h"

namespace MobiMesh {

/**
 * Concrete class for .tvpm64 mesh formats.
 */
class TrulyViewDependentPM64 : public TrulyViewDependentPM
{
public:

    TrulyViewDependentPM64();

    /**
     * Creates a TrulyViewDependentPM64 object from the given mesh filename
	 * @param[in] filename name of file containing tvpm64 data
	 * @throw FileReadErrorException if unable to read file
	 * @throw InvalidFileFormatException if file does not refer to a
     * TVPM64 data file
     */
    explicit TrulyViewDependentPM64(const std::string& filename);

    ~TrulyViewDependentPM64();

protected:
    // Note! We use 'unsigned int' (32-bit) or 'unsigned long (long)' (64-bit) as VertexID; (0.5.4)
    typedef unsigned long long      VertexID;
	
    typedef std::vector<VertexID>   VertexIDContainer; VertexIDContainer  vertex_vindex_;
	
    typedef struct ZsData {
        OpenMesh::Vec3f             p;
        OpenMesh::Vec3f             l_normal;
        OpenMesh::Vec3f             r_normal;
        unsigned int                fund_ccut_tree;
        unsigned int                fund_lcut_node;
        unsigned int                fund_rcut_node;
        int                         c_pos;
    } ZsData;
	
    unsigned int                    n_base_vertices_;
    unsigned int                    n_base_faces_;
    unsigned int                    n_details_;
	
    unsigned char                   tree_id_bits_;
    unsigned char                   tree_id_bitshift_; //                  64 - tree_id_bits_
    unsigned long long              tree_id_mask_;     // 0xFFFFFFFFFFFFFFFF >> tree_id_bits_
	
    // Note! We calculate the base offset and data size to access mapped files properly; (0.5.7)
    size_t                          s_base_, s_details_;

    ZsData loadVs (unsigned int pos);

    void get_active_cuts_fast (VPIter vpit  , unsigned int fund_lcut_tree, unsigned int fund_rcut_tree,
                               unsigned int fund_lcut_node, unsigned int fund_rcut_node,
                               int &vl, int &vr, int &v1vl, int &vrv1);

    void force_vsplit (VPIter vpit);
    void splt_index (int v0, int v1);
    void ecol_index (int v0, int v1);

    inline void swap_index (int v0, int v1) { std::swap(vertex_vindex_[v0], vertex_vindex_[v1]); }
    inline void resz_index (size_t vsize)   { vertex_vindex_.resize(vsize);                      }

    void calc_tree_id_bits (VertexID _n_roots) 
    {
        tree_id_bits_ = 0; while (_n_roots > ((VertexID) 1 << tree_id_bits_)) ++tree_id_bits_;
        tree_id_bitshift_ = sizeof(VertexID)*8   - tree_id_bits_;
        tree_id_mask_     =      ((VertexID)-1) >> tree_id_bits_;
    }

    void readVDPM(const char* _filename);
    
private:
    
    TrulyViewDependentPM64(const TrulyViewDependentPM64&);
    TrulyViewDependentPM64& operator=(const TrulyViewDependentPM64&);
};

}

#endif