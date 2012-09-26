//
//  TrulyViewDependentPM.h
//  MobiMesh
//
//  Created by Ong Yuh Shin on 5/18/11.
//  Copyright 2011 National University of Singapore. All rights reserved.
//

#ifndef __TRULY_VIEW_DEPENDENT_PM_H_
#define __TRULY_VIEW_DEPENDENT_PM_H_

#include <list>
#include <vector>
#include <map>

#include "IMesh.h"
#include "StatusInfo.h"

namespace MobiMesh {

/**
 * Abstract base class for a Truly View-dependent Progressive Mesh
 * For both .tvpm and .tvpm64 formats.
 */
class TrulyViewDependentPM : public IMesh
{
public:

    TrulyViewDependentPM();
    virtual ~TrulyViewDependentPM();

	///////////// Interface methods inherited from IMesh

    void load_mesh(const std::string& filename);

    const IMesh::PointContainer* points() const { return &points_; }

    const IMesh::VertexNormalContainer* vertex_normals() const { return NULL; }

    const IMesh::FaceContainer* faces() const { return NULL; }

    bool has_vertex_indices() const { return false; }

    const IMesh::PointContainer* face_vertices() const { return &vpoints_; }

    const IMesh::VertexNormalContainer* face_normals() const
    { return &vnormal_; }

    bool has_face_vertices() const { return true; }
    
    size_t num_faces() const { return faces_.size(); }

    ///////////// TrulyViewDependentPM specific methods

    /**
     * Adapts the mesh based on the vertex weights set.
     * This method has some preconditions. The edge collapse threshold need to
     * be strictly less than half of the vertex split threshold. This is
     * because each time an edge is collapsed, the remaining vertex's weight is
     * the sum of the weights of the 2 vertices forming the edge. If this sum
     * is greater than or equal to the vertex split threshold, the resulting
     * vertex needs to be split again. To avoid such ambiguity of whether to
     * split or collapse, we require that the edge collapse threshold be
     * strictly less than half the vertex split threshold.
     *
     * @param[in] max_refinements this method will perform no more than
     * max_refinements vertex splits and max_refinements edge collapses. Note
     * that this value is automatically reduced to MAX_ECVS_PER_FRAME by this
     * method if the user supplies a value above MAX_ECVS_PER_FRAME
     * @param[in] edge_collapse_threshold vertices whose weight are below or
     * equal to this threshold will be collapsed (when an edge is collapsed, 2
     * vertices are involved, the weight of both vertices need to be below this
     * threshold for the collapse to take place)
     * @param[in] vertex_split_threshold vertices whose weight are above or
     * equal to this threshold will be split
     *
     * @return true iff exactly max_refinements vertex splits or edge collapses
     * are performed, i.e. more refinements might be necessary based on the
     * supplied thresholds and current vertex weights.
     *
     * @throw InvalidThresholdException if
     * 2*edge_collapse_threshold >= vertex_split_threshold
     */
    bool adapt(int max_refinements=MAX_ECVS_PER_FRAME,
               unsigned int edge_collapse_threshold=1,
               unsigned int vertex_split_threshold=4);

    bool need_refined() const { return need_refined_; }

    /**
     * Adds the weight to the current weight of the vertex at the given index.
     * Use face_index_to_vertex_index to get the vertex_index for a particular
     * face. The assigned vertex weight is used to determine how the mesh
     * should be refined.
     */
    void add_vertex_weight(int vertex_index, unsigned int weight)
    { vertex_weight_[vertex_index] += weight; }
    void set_vertex_weight(int vertex_index, unsigned int weight)
    { vertex_weight_[vertex_index] = weight; }
    unsigned int get_vertex_weight(int vertex_index) const
    { return vertex_weight_[vertex_index]; }

    /**
     * Sets all the vertex weights to 0.
     */
    void reset_vertex_weights()
    {
        memset(&vertex_weight_[0], 0,
               vertex_weight_.size()*sizeof(unsigned int));
    }

    /**
     * A face i has 3 indices (3i, 3i+1, 3i+2) representing the 3 vertices of
     * the face. Each of these can be converted to the corresponding vertex
     * index to refer to the actual vertices of the face.
     * The vertex index can then be used to access/update the weight of the
     * particular vertex.
     *
     * @return vertex index of given face_index
     */
    int face_index_to_vertex_index(int face_index)
    { return draw_faces_[face_index]; }

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

    void unload_mesh();

    static const unsigned int     MAX_BUFF_FACES     = 200000;

protected:

    static const unsigned int     MAX_NEIGHBORS      = 300;
	static const unsigned int     MAX_ECVS_PER_FRAME = 2000;
	
    bool need_refined_;

    /// This flag is only set if mesh has been loaded successfully
    bool is_mesh_loaded_;

    /// This flag is only set if mmap() call was made and succeeded
    bool mmapped_;

    // Variables for mmap() and open()/pread().
    char* map_data; 
    long map_size;
    int fd;

    // Note! We cache some vertex neighborhood traverse results here; (0.6.1)
    int neighbors[MAX_NEIGHBORS];
    size_t ni;

    // *** MobiMesh equivalent: ViewDependentProgressiveMesh::typedef struct VParam/vplist_/vit_
    // Note! We keep view-dependent parameters of current vertex hierarchy (instead of explicit vertex front); (0.1.0)
    typedef struct ZParam {
        int vh; OpenMesh::Vec3f normal; unsigned int pos;
        std::list<ZParam>::iterator parent; 
        std::list<ZParam>::iterator lchild, rchild;
    } ZParam;
    typedef std::list<ZParam>           VPList; VPList vplist; 
    typedef std::list<ZParam>::iterator VPIter;

    // Note! We add 'eccount' and 'vscount' to limit the number of ec/vs every frame to improve the framerate; (0.5.7)
    // Note! We mark and re-allocate the free memory from edge-collapses to vertex-splits; (0.5.9)
    // Note! We unify the maximum number of ec/vs to maximize the efficiency of our memory re-allocation strategy; (0.6.2)
    // Note! We initialize the array size of free memory pointers according to the maximum number of ec/vs; (0.6.2)
    // Note! We know that every edge collapse create exactly 1 free vertex, 2 free faces, 3 free edges, 2 free vparams; (0.6.2)
    // Note! We need to know the possible overflow of vs counts due to force_vertex_split (neighbors) in every vs; (0.6.2)

    std::vector<size_t>   free_vert_;    int freevi;
    std::vector<size_t>   free_face_;    int freefi;
    std::vector<size_t>   free_edge_;    int freeei;
    std::vector<VPIter>   free_vplist_;  int freepi;

    // Note! We use OpenMesh geometry & connectivity data structure (not GPU-friendly); (0.1.0)
    typedef struct Vertex
	{
        int      halfedge_handle_;
        Vertex():halfedge_handle_(-1){}
    } Vertex;
    typedef struct Halfedge
	{
        int           face_handle_;
        int         vertex_handle_;
        int  next_halfedge_handle_;
        int  prev_halfedge_handle_; 
        Halfedge():face_handle_(-1),vertex_handle_(-1),next_halfedge_handle_(-1),prev_halfedge_handle_(-1){}
    } Halfedge;
    typedef struct Edge
	{
        Halfedge halfedges_[2];
    } Edge;
    typedef struct Face
	{
        int      halfedge_handle_;
        Face():halfedge_handle_(-1){}
    } Face;

    // Note! We use OpenMesh Vector to simulate the coordinate type and normal type; (0.1.0)
    // Note! We remark vertex neighbors while running a faster version of 'is_collapse_ok()'; (0.5.7)

    typedef OpenMesh::Vec3f  Point;
    typedef OpenMesh::Vec3f  Normal;

    typedef std::vector<Vertex>            VertexContainer;   VertexContainer   vertices_;
    typedef std::vector<Halfedge>          HfedgeContainer;   HfedgeContainer   hfedges_;
    typedef std::vector<Face>              FaceContainer;     FaceContainer     faces_;

    typedef std::vector<Point>             PointContainer;    PointContainer    points_;
    typedef std::vector<VPIter>            ZParamContainer;   ZParamContainer   vertex_params_;
    typedef std::vector<StatusInfo>        VStatusContainer;  VStatusContainer  vertex_status_;

    typedef std::vector<StatusInfo>        EStatusContainer;  EStatusContainer  edge_status_;
    typedef std::vector<StatusInfo>        FStatusContainer;  FStatusContainer  face_status_;

    // vertex_weight_[vertex_handle] gives the weight of the vertex identified by vertex_handle
    std::vector<unsigned int> vertex_weight_;

    // Face array, indexing the vertex of the face in the vertex container.
    std::vector<unsigned int>       draw_faces_;

    // Points and normals, vpoints_[3i+k] and vnormal_[3i+k] contain the point and normal
    // of vectors of the vertex at face i, which is made up of the vertices at indices
    // (3i+0,3i+1,3i+2)
    std::vector<OpenMesh::Vec3f>    vpoints_,vnormal_;

    inline size_t n_faces() const     { return      faces_.size();    }
    inline size_t n_vertices() const  { return   vertices_.size();    }
    inline size_t n_halfedges() const { return    hfedges_.size();    }
    inline size_t n_edges() const     { return    hfedges_.size()>>1; }

    inline bool vertices_empty() const  { return vertices_.empty(); }
    inline bool halfedges_empty() const { return  hfedges_.empty(); }
    inline bool edges_empty() const     { return  hfedges_.empty(); }
    inline bool faces_empty() const     { return    faces_.empty(); }

    int add_vertex(const Point& _p);
    int new_edge(int _start_vh, int _end_vh);
    int new_face();

    // Is halfedge _heh a boundary halfedge (is its face handle invalid) ?
    bool is_boundary_heh(int _heh) { return hfedges_[_heh].face_handle_ == -1; }
	
    // Is vertex _vh a boundary vertex ?
    bool is_boundary_vh (int _vh)
    {
        int heh = vertices_[_vh].halfedge_handle_; if (heh == -1) return true;
        int  fh = hfedges_[heh].face_handle_;  return fh == -1;
    }

	int face_handle(int _heh) { return hfedges_[_heh].face_handle_; }
    void set_face_handle(int _heh, int _fh) { hfedges_[_heh].face_handle_ = _fh; }

	int halfedge_handle_vh(int _vh) { return vertices_[_vh].halfedge_handle_; }
    void set_halfedge_handle_vh(int _vh, int _heh) { vertices_[_vh].halfedge_handle_ = _heh; }

	int halfedge_handle_fh(int _fh) { return faces_[_fh].halfedge_handle_; }
    void set_halfedge_handle_fh(int _fh, int _heh) { faces_[_fh].halfedge_handle_ = _heh; }

	int next_halfedge_handle(int _heh) { return hfedges_[ _heh].next_halfedge_handle_; }
	int prev_halfedge_handle(int _heh) { return hfedges_[ _heh].prev_halfedge_handle_; }
    void  set_next_halfedge_handle(int _heh, int _nheh)
    {
        hfedges_[ _heh].next_halfedge_handle_ = _nheh; 
		hfedges_[_nheh].prev_halfedge_handle_ =  _heh;
    }

    inline int  opposite_halfedge_handle (int _heh)      { return (_heh & 1) ? _heh-1 : _heh+1; }
    inline int  opp                      (int _heh)      { return (_heh & 1) ? _heh-1 : _heh+1; }

    int find_halfedge(int _start_vh, int _end_vh);
    void adjust_outgoing_halfedge(int _vh);

    template <class _Handle>
    struct NextCacheEntryT : public std::pair<_Handle, _Handle>{
	    typedef std::pair<_Handle, _Handle> Base;
	    NextCacheEntryT(_Handle _heh0, _Handle _heh1) : Base(_heh0, _heh1) { assert(_heh0 != -1); assert(_heh1 != -1); }
    };

    int add_face(unsigned int i0, unsigned int i1, unsigned int i2);

    inline static bool idxcomp (int i, int j) {return (i>j);} // decreasing order, e.g. 50,40,30,20,10..

    void garbage_collection_fast(bool _v, bool _e, bool _f);

    void vsplit_fast (VPIter vpit, const OpenMesh::Vec3f& p,
                      const OpenMesh::Vec3f& l_normal,
                      const OpenMesh::Vec3f& r_normal,
                      unsigned int l_pos, unsigned int r_pos,
                      int vl, int vr, int v1vl, int vrv1);

    bool ecol_fast (VPIter vpit);

    bool adaptive_refinement(int max_refinements,
                             unsigned int edge_collapse_threshold,
                             unsigned int vertex_split_threshold);

    virtual void readVDPM(const char * _filename) = 0;
    virtual void force_vsplit(VPIter vpit) = 0;
    virtual void splt_index (int v0, int v1) = 0;
    virtual void ecol_index (int v0, int v1) = 0;
    virtual void swap_index (int v0, int v1) = 0;
    virtual void resz_index (size_t   vsize) = 0;

private:

    TrulyViewDependentPM(const TrulyViewDependentPM&);
    TrulyViewDependentPM& operator=(const TrulyViewDependentPM&);

    bool is_active_;
    std::string mesh_filename_;
};

}

#endif