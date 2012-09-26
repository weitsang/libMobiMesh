#include <iostream>
#include <fstream>
#include <vector>
#include <list>
#include <map>
#include <algorithm>
#include <limits>
#include <exception>
#include <cmath>
#include <tr1/unordered_map>
#include <getopt.h>

#include "inc/VectorT.hh"
#include "inc/Timer.hh"

using namespace OpenMesh;

// Status bits used by the Status class. 
enum StatusBits {
    DELETED               = 1,    // Item has been deleted
    LOCKED                = 2,    // Item is locked.
    SELECTED              = 4,    // Item is selected.
    HIDDEN                = 8,    // Item is hidden.
    FEATURE               = 16,   // Item is a feature or belongs to a feature.
    TAGGED                = 32,   // Item is tagged.
    TAGGED2               = 64,   // Alternate bit for tagging an item.
    FIXEDNONMANIFOLD      = 128,  // Item was non-two-manifold and had to be fixed
    UNUSED                = 256   // Unused
};

// Note! We use 'unsigned char' instead of 'unsigned int' to save some memory; (0.5.6)
class StatusInfo {
public:  typedef unsigned char value_type; 
private: value_type status_;
public: 
    StatusInfo() : status_(0) {}
    unsigned char bits() const                    { return status_; }
    void      set_bits(unsigned char _bits)       { status_ = _bits; }
    bool    is_bit_set(unsigned char _s) const    { return (status_ & _s) > 0; }
    void       set_bit(unsigned char _s)          { status_ |= _s; }
    void     unset_bit(unsigned char _s)          { status_ &= ~_s; }
    void    change_bit(unsigned char _s, bool _b) { if (_b) status_ |= _s; else status_ &= ~_s; }

    bool     deleted()  const           { return is_bit_set(DELETED); }
    void set_deleted (bool _b)          { change_bit(DELETED, _b); }
    bool     locked()   const           { return is_bit_set(LOCKED); }
    void set_locked  (bool _b)          { change_bit(LOCKED, _b); }
    bool     selected() const           { return is_bit_set(SELECTED); }
    void set_selected(bool _b)          { change_bit(SELECTED, _b); }
    bool     hidden()   const           { return is_bit_set(HIDDEN); }
    void set_hidden  (bool _b)          { change_bit(HIDDEN, _b); }
    bool     feature()  const           { return is_bit_set(FEATURE); }
    void set_feature (bool _b)          { change_bit(FEATURE, _b); }
    bool     tagged()  const            { return is_bit_set(TAGGED); }
    void set_tagged  (bool _b)          { change_bit(TAGGED, _b); }
    bool     tagged2() const            { return is_bit_set(TAGGED2); }
    void set_tagged2 (bool _b)          { change_bit(TAGGED2, _b); }
    bool     fixed_nonmanifold() const  { return is_bit_set(FIXEDNONMANIFOLD); }
    void set_fixed_nonmanifold(bool _b) { change_bit(FIXEDNONMANIFOLD, _b); }
};


//////////   Out-of-Core View-Dependent Progressive Mesh   //////////

typedef unsigned int            VertexID;

unsigned int                    n_base_vertices_;
unsigned int                    n_base_faces_;
unsigned int                    n_details_;

unsigned char                   tree_id_bits_;
unsigned char                   tree_id_bitshift_; //          32 - tree_id_bits_
unsigned int                    tree_id_mask_;     // 0xFFFFFFFF >> tree_id_bits_

void calc_tree_id_bits(unsigned long long _n_roots)
{
    tree_id_bits_ = 0; while (_n_roots > ((VertexID) 1 << tree_id_bits_)) ++tree_id_bits_;
    tree_id_bitshift_ = sizeof(VertexID)*8   - tree_id_bits_;
    tree_id_mask_     =      ((VertexID)-1) >> tree_id_bits_;
}

// Note! We implement GPU-friendly face array and embed its updates to improve framerate; (0.5.8)
std::vector<unsigned int>       draw_faces_;

// Note! We keep view-dependent parameters of current vertex hierarchy (instead of explicit vertex front); (0.1.0)
// Note! We replace Vertex Hierarchy with Common Linked List (VParam Hierarchy); (0.1.0)
typedef struct VParam {
    int vh; VertexID vid, fund_lcut, fund_rcut;
    Vec3f normal; float radius, sin_square, mue_square, sigma_square;
    std::list<VParam>::iterator parent,lchild, rchild;
    VParam() : vh (-1), vid (0), fund_lcut (0), fund_rcut(0)
                , radius(0.0f), sin_square(0.0f), mue_square(0.0f), sigma_square(0.0f) {}
} VParam;
typedef std::list<VParam>               VPList; VPList  vplist;
typedef std::list<VParam>::iterator     VPIter;

//////////   Out-of-Core View-Dependent Progressive Mesh   //////////

// Note! We use OpenMesh geometry & connectivity data structure (not GPU-friendly); (0.1.0)
typedef struct Vertex{
	int      halfedge_handle_;
	Vertex():halfedge_handle_(-1){}
} Vertex;
typedef struct Halfedge{
	int           face_handle_;
	int         vertex_handle_;
	int  next_halfedge_handle_;
	int  prev_halfedge_handle_; 
	Halfedge():face_handle_(-1),vertex_handle_(-1),next_halfedge_handle_(-1),prev_halfedge_handle_(-1){}
} Halfedge;
typedef struct Edge{
	Halfedge halfedges_[2];
} Edge;
typedef struct Face{
	int  halfedge_handle_;
	Face():halfedge_handle_(-1){}
} Face;

// Note! We use OpenMesh Vector to simulate the coordinate type and normal type; (0.1.0)

typedef OpenMesh::Vec3f     Point;
typedef OpenMesh::Vec3f     Normal;

typedef struct VertInfo {
    VPIter  vparams_ptr;
} VertInfo;

typedef struct EdgeInfo {
    VPIter  vparams_ptr[2];
} EdgeInfo;

typedef std::vector<Vertex>                           VertexContainer;   VertexContainer       vertices_;
typedef std::vector<Edge>                               EdgeContainer;     EdgeContainer          edges_;
typedef std::vector<Face>                               FaceContainer;     FaceContainer          faces_;

typedef std::vector<Point>                             PointContainer;    PointContainer         points_;
typedef std::vector<Normal>                          VNormalContainer;  VNormalContainer  vertex_normal_;
typedef std::vector<StatusInfo>                      VStatusContainer;  VStatusContainer  vertex_status_;
typedef std::vector<VertInfo>                          VInfoContainer;    VInfoContainer  vertex_pminfo_;

typedef std::vector<StatusInfo>                      EStatusContainer;  EStatusContainer    edge_status_;
typedef std::vector<EdgeInfo>                          EInfoContainer;   EInfoContainer     edge_pminfo_;

typedef std::vector<Normal>                          FNormalContainer;  FNormalContainer    face_normal_;
typedef std::vector<StatusInfo>                      FStatusContainer;  FStatusContainer    face_status_;

// Note! We inline these utility functions to accelerate some queries; (0.5.7)

inline size_t    n_vertices()   { return vertices_.size(); }
inline size_t   n_halfedges()   { return  2*edges_.size(); }
inline size_t       n_edges()   { return    edges_.size(); }
inline size_t       n_faces()   { return    faces_.size(); }

inline bool  vertices_empty()   { return vertices_.empty(); }
inline bool halfedges_empty()   { return    edges_.empty(); }
inline bool     edges_empty()   { return    edges_.empty(); }
inline bool     faces_empty()   { return    faces_.empty(); }

// Note! We re-allocate the free vertex to vertex-splits before request a new one; (0.5.9)
int add_vertex(const Point& _p) 
{
	vertices_.push_back(Vertex());
	  points_.push_back(_p); 
	size_t nV = n_vertices();
	vertex_status_.resize(nV);
	vertex_normal_.resize(nV);
	vertex_pminfo_.resize(nV);
	return nV-1;
}
int new_edge(int _start_vh, int _end_vh)
{
	edges_.push_back(Edge());
	edge_status_.resize(n_edges());
	edge_pminfo_.resize(n_edges());
	int eh = n_edges()-1;
	edges_[eh].halfedges_[0].vertex_handle_ =   _end_vh;
	edges_[eh].halfedges_[1].vertex_handle_ = _start_vh;
	return eh<<1;
}
int new_face()
{
	faces_.push_back(Face());
	face_normal_.resize(n_faces());
	face_status_.resize(n_faces());
	return n_faces()-1;
}

 int to_vertex_handle          (int _heh)            { return edges_[_heh>>1].halfedges_[_heh&1].vertex_handle_;}

 int      face_handle          (int _heh)            { return edges_[_heh>>1].halfedges_[_heh&1].face_handle_; }
void  set_face_handle          (int _heh, int _fh)   {        edges_[_heh>>1].halfedges_[_heh&1].face_handle_ = _fh; }

 int      halfedge_handle_vh   (int _vh)             { return vertices_[_vh].halfedge_handle_; }
void  set_halfedge_handle_vh   (int _vh, int _heh)   {        vertices_[_vh].halfedge_handle_ = _heh; }

 int      halfedge_handle_fh   (int _fh)             { return faces_[_fh].halfedge_handle_; }
void  set_halfedge_handle_fh   (int _fh, int _heh)   {        faces_[_fh].halfedge_handle_ = _heh; }

 int      next_halfedge_handle (int _heh)            { return edges_[ _heh>>1].halfedges_[ _heh&1].next_halfedge_handle_; }
 int      prev_halfedge_handle (int _heh)            { return edges_[ _heh>>1].halfedges_[ _heh&1].prev_halfedge_handle_; }
void  set_next_halfedge_handle (int _heh, int _nheh) {        edges_[ _heh>>1].halfedges_[ _heh&1].next_halfedge_handle_ = _nheh; 
															  edges_[_nheh>>1].halfedges_[_nheh&1].prev_halfedge_handle_ =  _heh; }

 int  opposite_halfedge_handle (int _heh)            { return (_heh & 1) ? _heh-1 : _heh+1; }

// Is halfedge _heh a boundary halfedge (is its face handle invalid) ?
bool is_boundary_heh(int _heh) { return edges_[_heh>>1].halfedges_[_heh&1].face_handle_ == -1; }

// Is vertex _vh a boundary vertex ?
bool is_boundary_vh (int _vh)
{
	int heh = vertices_[_vh].halfedge_handle_; if (heh == -1) return true;
	int  fh = edges_[heh>>1].halfedges_[heh&1].face_handle_;  return fh == -1;
}

unsigned int valence_vh (int _vh)
{
    unsigned int count(0);
   
    int shh = vertices_[_vh].halfedge_handle_;
    int  hh = shh;
    do {if (hh == -1) break;
	    hh = edges_[hh >> 1].halfedges_[(hh&1) ^ 1].next_halfedge_handle_; ++count;
    } while (hh != shh);

    return count;
}

int find_halfedge(int _start_vh, int _end_vh)
{
	assert(_start_vh != -1 && _end_vh != -1);

	int shh = vertices_[_start_vh].halfedge_handle_;
	int  hh = shh;
	do {if (hh == -1) break;
		if ( edges_[hh >> 1].halfedges_[hh & 1].vertex_handle_ == _end_vh ) return hh;
		hh = edges_[hh >> 1].halfedges_[(hh&1) ^ 1].next_halfedge_handle_;
	} while (hh != shh);

	return -1;
}

void adjust_outgoing_halfedge(int _vh)
{
	int shh = vertices_[_vh].halfedge_handle_;
	int hh = shh;
	do {
		if (edges_[hh >> 1].halfedges_[hh & 1].face_handle_ == -1) {
			vertices_[_vh].halfedge_handle_ = hh;
			break; }
		hh = edges_[hh >> 1].halfedges_[(hh&1) ^ 1].next_halfedge_handle_;
	} while (hh != shh);
}

template <class _Handle>
struct NextCacheEntryT : public std::pair<_Handle, _Handle>{
	typedef std::pair<_Handle, _Handle> Base;
	NextCacheEntryT(_Handle _heh0, _Handle _heh1) : Base(_heh0, _heh1) { assert(_heh0 != -1); assert(_heh1 != -1); }
};

int add_face(unsigned int i0, unsigned int i1, unsigned int i2)
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
		if ( !is_boundary_vh(_vertex_handles[i]) )               { std::cerr << "PolyMeshT::add_face: complex vertex\n"; return -1; }
		halfedge_handles[i] = find_halfedge(_vertex_handles[i], _vertex_handles[ii]);
		is_new[i] = (halfedge_handles[i] == -1);
		if (!is_new[i] && !is_boundary_heh(halfedge_handles[i])) { std::cerr << "PolyMeshT::add_face: complex edge\n";   return -1; }
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
				if (boundary_next == inner_next)                 { std::cerr << "PolyMeshT::add_face: patch re-linking failed\n"; return -1; }
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

int vertex_split_fast(int v0, int v1, int vl, int vr)
{
    int v1vl = -1, vlv1 = -1, vrv1 = -1, v1vr = -1, v0v1 = -1;

    int sheh = vertices_[v1].halfedge_handle_;
    int  heh = sheh;
    do {
        int hehi = heh; if (hehi == -1) {v1vl=-1; vrv1=-1; break;}
        Halfedge& hehe  = edges_[hehi >> 1].halfedges_[hehi & 1];
        Halfedge& ohehe = edges_[hehi >> 1].halfedges_[(hehi&1) ^ 1];
        if (hehe.vertex_handle_ == vl) {v1vl=heh;}                           // v1vl = find_halfedge(v1, vl);
        if (hehe.vertex_handle_ == vr) {vrv1=opposite_halfedge_handle(heh);} // vrv1 = find_halfedge(vr, v1);
        heh = ohehe.next_halfedge_handle_;
    } while (heh != sheh);

    // build loop from halfedge v1->vl, vlv1 = insert_loop(v1vl);
    if (vl != -1)
    {
        // assert(v1vl != -1);

        // Note! We re-allocate the free edge to vertex-splits before request a new one; (0.5.9)
        int eidx = -1; edges_.push_back(Edge()); eidx = n_edges()-1;

        edges_[eidx].halfedges_[0].vertex_handle_ = v1;
        edges_[eidx].halfedges_[1].vertex_handle_ = vl;
        vlv1 = ((eidx<<1)+0);

        // Note! We re-allocate the free face to vertex-splits before request a new one; (0.5.9)
        int f1 = -1; faces_.push_back(Face()); draw_faces_.resize(n_faces()*3); f1 = n_faces()-1; 

        int            o0 = opposite_halfedge_handle(v1vl);
        int            o1 = opposite_halfedge_handle(vlv1);
        int            f0 = face_handle(v1vl);

        Halfedge&     he0 = edges_[v1vl >> 1].halfedges_[v1vl & 1];
        int           ph0 = he0.prev_halfedge_handle_;
        int           nh0 = he0.next_halfedge_handle_;
        Halfedge&    phe0 = edges_[ph0 >> 1].halfedges_[ph0 & 1];
        Halfedge&    nhe0 = edges_[nh0 >> 1].halfedges_[nh0 & 1];
        Halfedge&     he1 = edges_[vlv1 >> 1].halfedges_[vlv1 & 1];
        Halfedge&     oe0 = edges_[v1vl >> 1].halfedges_[(v1vl&1) ^ 1];
        Halfedge&     oe1 = edges_[vlv1 >> 1].halfedges_[(vlv1&1) ^ 1];

        phe0.next_halfedge_handle_ = o1; 
         oe1.prev_halfedge_handle_ = ph0;
         oe1.next_halfedge_handle_ = nh0; 
        nhe0.prev_halfedge_handle_ = o1;
         he1.next_halfedge_handle_ = v1vl;
         he0.prev_halfedge_handle_ = vlv1;
         he0.next_halfedge_handle_ = vlv1;
         he1.prev_halfedge_handle_ = v1vl;
        
        oe1.face_handle_ = f0;
        he0.face_handle_ = f1;
        he1.face_handle_ = f1;

                      faces_[f1].halfedge_handle_ = v1vl;
        if (f0 != -1) faces_[f0].halfedge_handle_ = o1;

        if(oe1.face_handle_ == -1) vertices_[v1].halfedge_handle_ = o1; // adjust_outgoing_halfedge of v1;
        if(oe0.face_handle_ == -1) vertices_[vl].halfedge_handle_ = o0; // adjust_outgoing_halfedge of vl;
    }

    // build loop from halfedge vr->v1, insert_loop(vrv1);
    if (vr != -1)
    {
        // assert(vrv1 != -1);

        // Note! We re-allocate the free edge to vertex-splits before request a new one; (0.5.9)
        int eidx = -1; edges_.push_back(Edge()); eidx = n_edges()-1;

        edges_[eidx].halfedges_[0].vertex_handle_ = vr;
        edges_[eidx].halfedges_[1].vertex_handle_ = v1;
        v1vr = ((eidx<<1)+0);
        
        // Note! We re-allocate the free face to vertex-splits before request a new one; (0.5.9)
        int f1 = -1; faces_.push_back(Face()); draw_faces_.resize(n_faces()*3); f1 = n_faces()-1; 

        int            o0 = opposite_halfedge_handle(vrv1);
        int            o1 = opposite_halfedge_handle(v1vr);
        int            f0 = face_handle(vrv1);

        Halfedge&     he0 = edges_[vrv1 >> 1].halfedges_[vrv1 & 1];
        int           ph0 = he0.prev_halfedge_handle_;
        int           nh0 = he0.next_halfedge_handle_;
        Halfedge&    phe0 = edges_[ph0 >> 1].halfedges_[ph0 & 1];
        Halfedge&    nhe0 = edges_[nh0 >> 1].halfedges_[nh0 & 1];
        Halfedge&     he1 = edges_[v1vr >> 1].halfedges_[v1vr & 1];
        Halfedge&     oe0 = edges_[vrv1 >> 1].halfedges_[(vrv1&1) ^ 1];
        Halfedge&     oe1 = edges_[v1vr >> 1].halfedges_[(v1vr&1) ^ 1];

        phe0.next_halfedge_handle_ = o1;  
         oe1.prev_halfedge_handle_ = ph0;
         oe1.next_halfedge_handle_ = nh0; 
        nhe0.prev_halfedge_handle_ = o1;
         he1.next_halfedge_handle_ = vrv1;
         he0.prev_halfedge_handle_ = v1vr;
         he0.next_halfedge_handle_ = v1vr;
         he1.prev_halfedge_handle_ = vrv1;
          
        oe1.face_handle_ = f0;
        he0.face_handle_ = f1;
        he1.face_handle_ = f1;

                      faces_[f1].halfedge_handle_ = vrv1;
        if (f0 != -1) faces_[f0].halfedge_handle_ = o1;
        
        if(oe1.face_handle_ == -1) vertices_[vr].halfedge_handle_ = o1; // adjust_outgoing_halfedge of vr;
        if(oe0.face_handle_ == -1) vertices_[v1].halfedge_handle_ = o0; // adjust_outgoing_halfedge of v1;
    }

    // handle boundary cases
    if (vl == -1) vlv1 = prev_halfedge_handle(halfedge_handle_vh(v1));
    if (vr == -1) vrv1 = prev_halfedge_handle(halfedge_handle_vh(v1));
    // assert(vlv1 != -1 && vrv1 != -1 && v1 == to_vertex_handle(vrv1));

    // split vertex v1 into edge v0v1, v0v1 = insert_edge(v0, vlv1, vrv1);
    
    // Note! We re-allocate the free edge to vertex-splits before request a new one; (0.5.9)
    int eidx = -1; edges_.push_back(Edge()); eidx = n_edges()-1; 

    edges_[eidx].halfedges_[0].vertex_handle_ = v1;
    edges_[eidx].halfedges_[1].vertex_handle_ = v0;
    v0v1 = ((eidx<<1)+0);
    int v1v0 = opposite_halfedge_handle(v0v1);

    Halfedge&     he01 = edges_[v0v1 >> 1].halfedges_[v0v1 & 1];
    Halfedge&     he10 = edges_[v1v0 >> 1].halfedges_[v1v0 & 1];
    Halfedge&     hel1 = edges_[vlv1 >> 1].halfedges_[vlv1 & 1];
    Halfedge&     her1 = edges_[vrv1 >> 1].halfedges_[vrv1 & 1];
    int           nhl1 = hel1.next_halfedge_handle_;
    int           nhr1 = her1.next_halfedge_handle_;
    Halfedge&    nhel1 = edges_[nhl1 >> 1].halfedges_[nhl1 & 1];
    Halfedge&    nher1 = edges_[nhr1 >> 1].halfedges_[nhr1 & 1]; 

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

    if (he01.face_handle_ != -1) faces_[he01.face_handle_].halfedge_handle_ = v0v1;
    if (he10.face_handle_ != -1) faces_[he10.face_handle_].halfedge_handle_ = v1v0;

    // we do not always need to adjust the outgoing (boundary) halfedge of v1.
    Vertex& vv0 = vertices_[v0];
    Vertex& vv1 = vertices_[v1];
    int   hehv1 = vv1.halfedge_handle_; vv0.halfedge_handle_ = v0v1;
            heh = v0v1; int df = -1, pdf = he01.face_handle_;
    do {
        Halfedge&  hehe = edges_[heh >> 1].halfedges_[heh & 1];
        Halfedge& ohehe = edges_[heh >> 1].halfedges_[(heh&1) ^ 1];
        if (heh == hehv1)            {vv1.halfedge_handle_ = v1v0;} // adjust outgoing halfedge of v1
        if (hehe.face_handle_ == -1) {vv0.halfedge_handle_ = heh;}  // adjust outgoing halfedge of v0
        ohehe.vertex_handle_ = v0;

        // Note! We embed the updates of GPU-friendly 'draw_faces_' into vertex-splits; (0.5.8)
        // Note! We fix the segmentation fault while updating draw_faces_ on neptune_nonmanifold.vpm; (0.5.9)
        df = ohehe.face_handle_; int vd = hehe.vertex_handle_;
        if ( df >= 0) draw_faces_[ df*3  ] = v0; // ohehe.vertex_handle_;
        if (pdf >= 0) draw_faces_[pdf*3+1] = vd; //  hehe.vertex_handle_;
        if ( df >= 0) draw_faces_[ df*3+2] = vd; //  hehe.vertex_handle_;
        pdf = df;

        heh = ohehe.next_halfedge_handle_;
    } while (heh != v0v1);

    // Note! We move these to 'adaptive_refinement()' to speed up a single vertex-split; (0.1.0)
    // eprops_resize(n_edges());
    // hprops_resize(n_halfedges());
    // fprops_resize(n_faces());

    // make sure the outgoing halfedge adjustment conforms to the assumption (only one boundary per vertex)
    // for(VertexOHalfedgeIter voh_it(voh_iter(v0)); voh_it; ++voh_it) {
    //     if (is_boundary(voh_it.handle())) assert(voh_it.handle()==vertices_[v0].halfedge_handle_);}
    // for(VertexOHalfedgeIter voh_it(voh_iter(v1)); voh_it; ++voh_it) {
    //     if (is_boundary(voh_it.handle())) assert(voh_it.handle()==vertices_[v1].halfedge_handle_);}

    return v0v1;
}

void collapse_fast(int _hh)
{
	// remove edge, collapse_edge(v0v1)

	int  v0v1 = _hh;                            Halfedge&  he01 = edges_[v0v1 >> 1].halfedges_[v0v1 & 1];
	int  v1v0 = opposite_halfedge_handle(v0v1); Halfedge&  he10 = edges_[v1v0 >> 1].halfedges_[v1v0 & 1];

	int    hn = he01.next_halfedge_handle_;     Halfedge&   hne = edges_[  hn >> 1].halfedges_[  hn & 1];
	int    hp = he01.prev_halfedge_handle_;     Halfedge&   hpe = edges_[  hp >> 1].halfedges_[  hp & 1];
	int    on = he10.next_halfedge_handle_;     Halfedge&   one = edges_[  on >> 1].halfedges_[  on & 1];
	int    op = he10.prev_halfedge_handle_;     Halfedge&   ope = edges_[  op >> 1].halfedges_[  op & 1];
	
	int    fh = he01.face_handle_;
	int    fo = he10.face_handle_;
	int    v1 = he01.vertex_handle_;
	int    v0 = he10.vertex_handle_;

	// Note! Make sure the outgoing halfedge adjustment conforms to the assumption (only one boundary per vertex);
	// for(VertexOHalfedgeIter voh_it(voh_iter(v0)); voh_it; ++voh_it) {
	//     if (is_boundary(voh_it.handle())) assert(voh_it.handle()==vertices_[v0].halfedge_handle_);}
	// for(VertexOHalfedgeIter voh_it(voh_iter(v1)); voh_it; ++voh_it) {
	//     if (is_boundary(voh_it.handle())) assert(voh_it.handle()==vertices_[v1].halfedge_handle_);}

	// Note! We do not always need to adjust the outgoing (boundary) halfedge of v1;
	Vertex& vv0 = vertices_[v0]; Vertex& vv1 = vertices_[v1]; 
	int heh = v0v1; int df = -1, pdf = he01.face_handle_;
	do {
		Halfedge&  hehe = edges_[heh >> 1].halfedges_[heh & 1];
		Halfedge& ohehe = edges_[heh >> 1].halfedges_[(heh&1) ^ 1];
		ohehe.vertex_handle_ = v1;

		// Note! We embed the updates of GPU-friendly 'draw_faces_' into edge-collapses; (0.5.8)
		// Note! We fix the segmentation fault while updating draw_faces_ on neptune_nonmanifold.vpm; (0.5.9)
		df = ohehe.face_handle_; int vd = hehe.vertex_handle_;
		if ( df >= 0) draw_faces_[ df*3  ] = v1; // ohehe.vertex_handle_;
		if (pdf >= 0) draw_faces_[pdf*3+1] = vd; //  hehe.vertex_handle_;
		if ( df >= 0) draw_faces_[ df*3+2] = vd; //  hehe.vertex_handle_;
		pdf = df;

		heh = ohehe.next_halfedge_handle_;
	} while (heh != v0v1);

	// adjust outgoing halfedge of v1
	Halfedge& vhe0 = edges_[vv0.halfedge_handle_ >> 1].halfedges_[vv0.halfedge_handle_ & 1];
	if (vv0.halfedge_handle_ != v0v1 && vhe0.face_handle_ == -1) vv1.halfedge_handle_ = vv0.halfedge_handle_; 
	if (vv1.halfedge_handle_ == v1v0)                            vv1.halfedge_handle_ = hn;
	// assert(to_vertex_handle(opposite_halfedge_handle(vv1.halfedge_handle_)) == v1);

	hpe.next_halfedge_handle_ = hn;
	hne.prev_halfedge_handle_ = hp;
	ope.next_halfedge_handle_ = on;
	one.prev_halfedge_handle_ = op;

	if (fh != -1)  faces_[fh].halfedge_handle_ = hn;
	if (fo != -1)  faces_[fo].halfedge_handle_ = on;

	// Note! We mark the free edge and vertex from edge-collapses; (0.5.9)
	vertices_[v0].halfedge_handle_ = -1; // invalidate(), set_isolated(v0);
	  edge_status_[v0v1>>1].set_deleted(true);  
	vertex_status_[v0]     .set_deleted(true);

	// remove loops, collapse_loop(next_halfedge_handle(hn)) 
	if (hne.next_halfedge_handle_ == hp && hpe.prev_halfedge_handle_ == hn) {

		int   h0 = hp;                           Halfedge&   he0 = hpe;
		int   h1 = hn;                           Halfedge&   he1 = hne;
		
		int   o0 = opposite_halfedge_handle(h0); Halfedge&   oe0 = edges_[ o0 >> 1].halfedges_[ o0 & 1];
		int  op0 = oe0.prev_halfedge_handle_;    Halfedge&  ope0 = edges_[op0 >> 1].halfedges_[op0 & 1];
		int  on0 = oe0.next_halfedge_handle_;    Halfedge&  one0 = edges_[on0 >> 1].halfedges_[on0 & 1];
		
		int   v0 = he0.vertex_handle_; // v1
		int   v1 = he1.vertex_handle_; // vl
		int   fh = he0.face_handle_;
		int   fo = oe0.face_handle_;

		 he1.next_halfedge_handle_ = on0;
		one0.prev_halfedge_handle_ =  h1;
		 he1.prev_halfedge_handle_ = op0;
		ope0.next_halfedge_handle_ =  h1;

		he1.face_handle_ = fo;

		if (vertices_[v0].halfedge_handle_ == o0) vertices_[v0].halfedge_handle_ = h1;  // adjust outgoing halfedge of v0
		if (vertices_[v1].halfedge_handle_ == h0) vertices_[v1].halfedge_handle_ = on0; // adjust outgoing halfedge of v1

		if (fo != -1 && faces_[fo].halfedge_handle_ == o0) {faces_[fo].halfedge_handle_ = h1;}

		// Note! We mark the free face and edge from edge-collapses; (0.5.9)
		if (fh != -1)  {
			faces_[fh].halfedge_handle_ = -1;
			face_status_[fh].set_deleted(true); }
		 edge_status_[h0>>1].set_deleted(true); 
	}

	// remove loops, collapse_loop(on)
	if (one.next_halfedge_handle_ == op && ope.prev_halfedge_handle_ == on) {

		int   h0 = on;                           Halfedge&   he0 = one;
		int   h1 = op;                           Halfedge&   he1 = ope;
		
		int   o0 = opposite_halfedge_handle(h0); Halfedge&   oe0 = edges_[ o0 >> 1].halfedges_[ o0 & 1];
		int  op0 = oe0.prev_halfedge_handle_;    Halfedge&  ope0 = edges_[op0 >> 1].halfedges_[op0 & 1];
		int  on0 = oe0.next_halfedge_handle_;    Halfedge&  one0 = edges_[on0 >> 1].halfedges_[on0 & 1];

		int   v0 = he0.vertex_handle_; // vr
		int   v1 = he1.vertex_handle_; // v1
		int   fh = he0.face_handle_;
		int   fo = oe0.face_handle_;

		 he1.next_halfedge_handle_ = on0;
		one0.prev_halfedge_handle_ =  h1;
		 he1.prev_halfedge_handle_ = op0;
		ope0.next_halfedge_handle_ =  h1;

		he1.face_handle_ = fo;

		if (vertices_[v0].halfedge_handle_ == o0) vertices_[v0].halfedge_handle_ = h1;  // adjust outgoing halfedge of v0
		if (vertices_[v1].halfedge_handle_ == h0) vertices_[v1].halfedge_handle_ = on0; // adjust outgoing halfedge of v1

		if (fo != -1 && faces_[fo].halfedge_handle_ == o0) {faces_[fo].halfedge_handle_ = h1;}

		// Note! We mark the free face and edge from edge-collapses; (0.5.9)
		if (fh != -1)  {
			faces_[fh].halfedge_handle_ = -1;
			face_status_[fh].set_deleted(true); }
		 edge_status_[h0>>1].set_deleted(true); 
	}
}

void garbage_collection_fast (bool _v, bool _e, bool _f)
{
    int i0, i1, nV(n_vertices()), nE(n_edges()), nF(n_faces());

    // remove deleted vertices
    if (_v && n_vertices() > 0)
    {
        i0=0;  i1=nV-1;
        while (1) {
            // find 1st deleted and last un-deleted
            while (!vertex_status_[i0].deleted() && i0 < i1)  ++i0;
            while ( vertex_status_[i1].deleted() && i0 < i1)  --i1;
            if (i0 >= i1) break;

            // update handles of halfedges
            int  sheh = vertices_[i1].halfedge_handle_; Halfedge& shehe = edges_[sheh >> 1].halfedges_[sheh & 1];
            int   heh = sheh; int df = -1, pdf = shehe.face_handle_;
            do {
                Halfedge&  hehe = edges_[heh >> 1].halfedges_[heh & 1];
                Halfedge& ohehe = edges_[heh >> 1].halfedges_[(heh&1) ^ 1];
                ohehe.vertex_handle_ = i0;

                // Note! We embed the updates of GPU-friendly 'draw_faces_' into garbage_collection; (0.5.8)
                // Note! We fix the segmentation fault while updating draw_faces_ on neptune_nonmanifold.vpm; (0.5.9)
                df = ohehe.face_handle_; int vd = hehe.vertex_handle_;
                if ( df >= 0) draw_faces_[ df*3  ] = i0; // ohehe.vertex_handle_;
                if (pdf >= 0) draw_faces_[pdf*3+1] = vd; //  hehe.vertex_handle_;
                if ( df >= 0) draw_faces_[ df*3+2] = vd; //  hehe.vertex_handle_;
                pdf = df;

                heh = ohehe.next_halfedge_handle_;
            } while (heh != sheh);

            std::swap(     vertices_[i0],      vertices_[i1]);
            std::swap(       points_[i0],        points_[i1]);
            std::swap(vertex_normal_[i0], vertex_normal_[i1]);
            std::swap(vertex_status_[i0], vertex_status_[i1]);
            std::swap(vertex_pminfo_[i0], vertex_pminfo_[i1]);
        };
             vertices_.resize(vertex_status_[i0].deleted() ? i0 : i0+1);
               points_.resize(n_vertices());
        vertex_normal_.resize(n_vertices());
        vertex_status_.resize(n_vertices());
        vertex_pminfo_.resize(n_vertices());
    }
    // remove deleted faces
    if (_f && n_faces() > 0)
    {
        i0=0;  i1=nF-1;
        while (1)
        {
            // find 1st deleted and last un-deleted
            while (!face_status_[i0].deleted() && i0 < i1)  ++i0;
            while ( face_status_[i1].deleted() && i0 < i1)  --i1;
            if (i0 >= i1) break;

            // update handles of halfedges
            int   fh = faces_[i1].halfedge_handle_; Halfedge&  fhe = edges_[ fh >> 1].halfedges_[ fh & 1];  fhe.face_handle_ = i0;
            int  pfh =   fhe.prev_halfedge_handle_; Halfedge& pfhe = edges_[pfh >> 1].halfedges_[pfh & 1]; pfhe.face_handle_ = i0;
            int  nfh =   fhe.next_halfedge_handle_; Halfedge& nfhe = edges_[nfh >> 1].halfedges_[nfh & 1]; nfhe.face_handle_ = i0;
            
            std::swap(      faces_[i0],       faces_[i1]);
            std::swap(face_normal_[i0], face_normal_[i1]);
            std::swap(face_status_[i0], face_status_[i1]);

            // Note! We swap the elements of GPU-friendly 'draw_faces_' in garbage_collection; (0.5.8)
            draw_faces_[i0*3  ] = pfhe.vertex_handle_; // std::swap(drawfaces[i0*3  ],  drawfaces[i1*3  ]);
            draw_faces_[i0*3+1] =  fhe.vertex_handle_; // std::swap(drawfaces[i0*3+1],  drawfaces[i1*3+1]);
            draw_faces_[i0*3+2] = nfhe.vertex_handle_; // std::swap(drawfaces[i0*3+2],  drawfaces[i1*3+2]);
        };
              faces_.resize(face_status_[i0].deleted() ? i0 : i0+1);
        face_normal_.resize(n_faces());
        face_status_.resize(n_faces());

         draw_faces_.resize(n_faces()*3);
    }
    // remove deleted edges
    if (_e && n_edges() > 0)
    {
        i0=0; i1=nE-1;
        while (1)
        {
            // find 1st deleted and last un-deleted
            while (!edge_status_[i0].deleted() && i0 < i1)  ++i0;
            while ( edge_status_[i1].deleted() && i0 < i1)  --i1;
            if (i0 >= i1) break;

            int     h0 = 2*i0  ; Halfedge&  he0 = edges_[i1].halfedges_[0];
            int     h1 = 2*i0+1; Halfedge&  he1 = edges_[i1].halfedges_[1];
            // update handles of vertices
            Vertex& v0 = vertices_[he0.vertex_handle_];
            Vertex& v1 = vertices_[he1.vertex_handle_];
            if (v0.halfedge_handle_ == 2*i1+1) v0.halfedge_handle_ = h1;
            if (v1.halfedge_handle_ == 2*i1  ) v1.halfedge_handle_ = h0;
            // update handles of halfedges (next/prev)
            int    ph0 = he0.prev_halfedge_handle_;
            int    nh0 = he0.next_halfedge_handle_;
            edges_[ph0 >> 1].halfedges_[ph0 & 1].next_halfedge_handle_ = h0;
            edges_[nh0 >> 1].halfedges_[nh0 & 1].prev_halfedge_handle_ = h0;
            int    ph1 = he1.prev_halfedge_handle_;
            int    nh1 = he1.next_halfedge_handle_;
            edges_[ph1 >> 1].halfedges_[ph1 & 1].next_halfedge_handle_ = h1;
            edges_[nh1 >> 1].halfedges_[nh1 & 1].prev_halfedge_handle_ = h1;
            // update handles of faces
            Face&   f0 = faces_[he0.face_handle_];
            Face&   f1 = faces_[he1.face_handle_];
            if (f0.halfedge_handle_ == 2*i1  ) f0.halfedge_handle_ = h0;
            if (f1.halfedge_handle_ == 2*i1+1) f1.halfedge_handle_ = h1;
            
            std::swap(      edges_[i0],       edges_[i1]);
            std::swap(edge_status_[i0], edge_status_[i1]);
        };
              edges_.resize(edge_status_[i0].deleted() ? i0 : i0+1);
        edge_status_.resize(n_edges());
        edge_pminfo_.resize(n_edges());
    }
}

void update_normals(std::vector<unsigned int>& faces)
{
    int vertices_size = vertices_.size();
    std::vector<unsigned int> vertex_faces; vertex_faces.resize(vertices_size);
    for (int vi = 0; vi < vertices_size; ++vi) {
        vertex_normal_[vi] = Normal(0,0,0);
        vertex_faces  [vi] = 0;
    }
    int faces_size = faces_.size(); faces.resize(3*faces_size);
    for (int fi = 0; fi < faces_size; ++fi) {
        Face& f = faces_[fi];
        int heh = f.halfedge_handle_;
        Halfedge& he = edges_[heh>>1].halfedges_[heh&1];
        int nheh = he.next_halfedge_handle_;
        int pheh = he.prev_halfedge_handle_;
        Halfedge& nhe = edges_[nheh>>1].halfedges_[nheh&1];
        Halfedge& phe = edges_[pheh>>1].halfedges_[pheh&1];
        Point& p0 = points_[phe.vertex_handle_];
        Point& p1 = points_[ he.vertex_handle_];
        Point& p2 = points_[nhe.vertex_handle_]; 
        Normal p1p0 = p0 - p1;
        Normal p1p2 = p2 - p1;
        Normal n = p1p2 % p1p0; float norm = n.length();
        face_normal_[fi] = ( (norm != float(0)) ? ((n *= (float(1)/norm)),n) : Normal(0,0,0) );

        vertex_normal_[phe.vertex_handle_] += face_normal_[fi]; vertex_faces[phe.vertex_handle_] += 1;
        vertex_normal_[ he.vertex_handle_] += face_normal_[fi]; vertex_faces[ he.vertex_handle_] += 1;
        vertex_normal_[nhe.vertex_handle_] += face_normal_[fi]; vertex_faces[nhe.vertex_handle_] += 1;

        faces[3*fi  ] = phe.vertex_handle_;
        faces[3*fi+1] =  he.vertex_handle_;
        faces[3*fi+2] = nhe.vertex_handle_;
    }
    for (int vi = 0; vi < vertices_size; ++vi) {
        if (vertex_faces[vi] != 0) vertex_normal_[vi] /= vertex_faces[vi];
    }
}


typedef struct PMInfo { Vec3f p0; int v0, v1, vl, vr; } PMInfo;

typedef std::vector<PMInfo>                     PMInfoContainer;
typedef std::vector<PMInfo>::iterator           PMInfoIter;
typedef std::vector<VPIter>                     VPIterContainer;
typedef std::vector<Vec3f>                      ResidualContainer;

PMInfoContainer                                 pminfos_;
PMInfoIter                                      pmiter_;
unsigned int                                    n_current_res_;
unsigned int                                    n_max_res_;
bool                                            verbose = false;

void readPM                     (const char * _filename); // open progressive mesh
void lovePM                     ();
void savePM                     (const char * _filename); // save view-dependent progressive mesh

void refinePM(unsigned int _n);      //  refine mesh   up to _n vertices
void coarsePM(unsigned int _n);      // coarsen mesh down to _n vertices

void get_leaf_node_handles      (VPIter vp, VPIterContainer &leaf_list);
void compute_bounding_box       (VPIter vp, VPIterContainer &leaf_list);
void compute_cone_of_normals    (VPIter vp, VPIterContainer &leaf_list);
void compute_screen_space_error (VPIter vp, VPIterContainer &leaf_list);
void compute_mue_sigma          (VPIter vp, ResidualContainer &residuals);

// Vec3f point2triangle_residual   (const Vec3f &p, const Vec3f tri[3], float &s, float &t);
Vec3f point2triangle_residual(const Vec3f &p, const Vec3f& B, const Vec3f& E0, const Vec3f& E1, float &s, float &t);


void usage_and_exit(int xcode)
{
    using namespace std;
    cout << "Usage: vdpmanalyzer [-h] [-o output.vpm] input.pm\n";
    exit(xcode);
}

inline std::string& replace_extension( std::string& _s, const std::string& _e )
{
    std::string::size_type dot = _s.rfind(".");
    if (dot == std::string::npos) { _s += "." + _e            ; }
    else                          { _s = _s.substr(0,dot+1)+_e; }
    return _s;
}

inline std::string basename(const std::string& _f)
{
    std::string::size_type dot = _f.rfind("/");
    if (dot == std::string::npos) return _f;
    return _f.substr(dot+1, _f.length()-(dot+1));
}

int main (int argc, char **argv)
{
    int           c;
    std::string   ifname;
    std::string   ofname;

    while ( (c=getopt(argc, argv, "o:"))!=-1 ) {
        switch(c) {
        case 'v': verbose = true;   break;
        case 'o': ofname = optarg;  break;
        case 'h': usage_and_exit(0);
        default:  usage_and_exit(1);
    }}

    if (optind >= argc) usage_and_exit(1);

    ifname = argv[optind];

    if (ofname == "." || ofname == ".." )
    ofname += "/" + basename(ifname);
    std::string vpmfname = ofname.empty() ? ifname : ofname;
    replace_extension(vpmfname, "vpm");

    if ( ifname.empty() || vpmfname.empty() ) usage_and_exit(1);

    try {
        readPM (ifname.c_str());
        lovePM();
        savePM (vpmfname.c_str());
    }
    catch ( std::bad_alloc&   ) {std::cerr << "Error: out of memory!\n" << std::endl;   return 1;}
    catch ( std::exception& x ) {std::cerr << "Error: " << x.what() << std::endl;       return 1;}
    catch ( ... )               {std::cerr << "Fatal! Unknown error!\n";                return 1;}
    return 0;
}

VPIterContainer rootvec;

void readPM (const char * _filename)
{

    char                        fileformat[10];
    Vec3f                       p;
    int                         v1,vl,vr,fvi[3];

    int                         vertex_handle;
    VertexID                    node_index;

    std::ifstream ifs(_filename, std::ios::binary);

    ifs.read(fileformat, 8); fileformat[8]='\0';
    if (std::string(fileformat) != std::string("ProgMesh")) 
    { std::cerr << "Wrong file format.\n"; ifs.close(); exit(1); }

    ifs.read( (char*)&n_base_vertices_, sizeof(n_base_vertices_) );
    ifs.read( (char*)&n_base_faces_   , sizeof(n_base_faces_) );
    ifs.read( (char*)&n_details_      , sizeof(n_details_) );

    calc_tree_id_bits(n_base_vertices_);

    for (size_t i=0; i<n_base_vertices_; ++i) 
    {
        ifs.read( (char*)&p , sizeof(p) ); 

        vertex_handle   = add_vertex(p);
        node_index      = ( (VertexID) i << tree_id_bitshift_) | 1;

        VParam vp;
        vp.vh           = vertex_handle; assert(vertex_handle == (int) i);
        vp.vid          = node_index;
        vp.parent       = vplist.end(); 
        vp.lchild       = vplist.end(); 
        vp.rchild       = vplist.end();

        vertex_pminfo_[vertex_handle].vparams_ptr = vplist.insert(vplist.end(),vp);

        assert(vertex_handle == (int) i);
        rootvec.push_back(vertex_pminfo_[i].vparams_ptr);

    }
    for (size_t i=0; i<n_base_faces_; ++i)
    {
        ifs.read( (char*)&fvi[0] , sizeof(fvi[0]) );
        ifs.read( (char*)&fvi[1] , sizeof(fvi[1]) );
        ifs.read( (char*)&fvi[2] , sizeof(fvi[2]) );

        add_face(fvi[0],fvi[1],fvi[2]);
    }
    for (size_t i=0; i<n_details_; ++i)
    {
        ifs.read( (char*)&p  , sizeof(p)  );
        ifs.read( (char*)&v1 , sizeof(v1) ); 
        ifs.read( (char*)&vl , sizeof(vl) ); 
        ifs.read( (char*)&vr , sizeof(vr) ); 

        PMInfo pminfo;
        pminfo.p0 = p;
        pminfo.v0 = add_vertex(p);
        pminfo.v1 = v1;
        pminfo.vl = vl;
        pminfo.vr = vr;
        pminfos_.push_back(pminfo);

        VPIter            pvp = vertex_pminfo_[pminfo.v1].vparams_ptr;

        VertexID    parent_id = pvp->vid;
        VertexID      tree_id = parent_id >> tree_id_bitshift_;
        VertexID      node_id = parent_id  & tree_id_mask_    ;
        VertexID    lchild_id = ( tree_id << tree_id_bitshift_) | (2*node_id  );
        VertexID    rchild_id = ( tree_id << tree_id_bitshift_) | (2*node_id+1);

        VParam lvp, rvp;
        lvp.vh = pminfo.v0; lvp.vid = lchild_id; lvp.parent = pvp; lvp.lchild = vplist.end(); lvp.rchild = vplist.end();
        rvp.vh = pminfo.v1; rvp.vid = rchild_id; rvp.parent = pvp; rvp.lchild = vplist.end(); rvp.rchild = vplist.end();

        VPIter lchild = vplist.insert(vplist.end(),lvp); pvp->lchild = lchild; vertex_pminfo_[pminfo.v0].vparams_ptr = lchild;
        VPIter rchild = vplist.insert(vplist.end(),rvp); pvp->rchild = rchild; vertex_pminfo_[pminfo.v1].vparams_ptr = rchild;
    }
    ifs.close();

    assert(rootvec.size() == n_base_vertices_);

    // recover mapping between basemesh vertices to roots of vertex hierarchy
    for (size_t i = 0; i < n_base_vertices_; ++i)
    {
        vertex_pminfo_[i].vparams_ptr = rootvec[i];
    }

    pmiter_ = pminfos_.begin();
    n_current_res_ = 0;
    n_max_res_ = n_details_;

    // update face and vertex normals
    update_normals(draw_faces_);

    // bounding box
    Vec3f bbMin = points_[0], bbMax = points_[0];

    for (size_t i = 0; i < n_vertices(); i++) { bbMin.minimize(points_[i]); bbMax.maximize(points_[i]); }

    std::cerr << n_vertices() << " vertices, "
              << n_edges()    << " edge, "
              << n_faces()    << " faces, "
              << n_details_   << " detail vertices\n";
}

// Note! We load view-dependent parameters from secondary storage when necessary;
typedef struct VsData32 {
    Vec3f                       p;
    unsigned int                node_index;
    unsigned int                fund_lcut_index; //  left cut neighbor
    unsigned int                fund_rcut_index; // right cut neighbor
    float                       l_radius;
    Vec3f                       l_normal;
    float                       l_sin_square;
    float                       l_mue_square;
    float                       l_sigma_square;
    float                       r_radius;
    Vec3f                       r_normal;
    float                       r_sin_square;
    float                       r_mue_square;
    float                       r_sigma_square;
    unsigned int                l_pos;
    unsigned int                r_pos;
} VsData32;

std::tr1::unordered_map<VertexID, unsigned int> mobivmem;

void statPM()
{
    // write base mesh  
    coarsePM(0);
    garbage_collection_fast ( false, true, true );

    // write progressive detail (vertex hierarchy)
    for (size_t i=0; i<n_details_; ++i)
    {
        PMInfo pminfo = *pmiter_;
        VPIter vp = vertex_pminfo_[pminfo.v1].vparams_ptr;
        mobivmem[vp->vid] = i+1;
        refinePM(i);
    }
}

void savePM (const char * _filename)
{    
    unsigned int                            fvi[3], pos;
    Vec3f                                   p;
    Vec3f                                   normal;
    float                                   radius, sin_square, mue_square, sigma_square;
    int                                     vh;

    std::map<int, unsigned int>             handle2index_map;

    statPM();

    std::ofstream ofs(_filename, std::ios::binary);
    if (!ofs) {std::cerr << "write error\n"; exit(1);}
    
    ofs << "VDProgMesh"; // write header

    ofs.write( (char*)&n_base_vertices_, sizeof(n_base_vertices_) );
    ofs.write( (char*)&n_base_faces_   , sizeof(n_base_faces_   ) );
    ofs.write( (char*)&n_details_      , sizeof(n_details_      ) );

    // write base mesh  
    coarsePM(0);
    garbage_collection_fast( false, true, true );

    assert(rootvec.size()==n_base_vertices_);
    for (size_t i = 0; i < n_base_vertices_; ++i)
    {
        VPIter vp = rootvec[i]; 

        vh            = vp->vh;

        p             = points_[vh];
        radius        = vp->radius;
        normal        = vp->normal;
        sin_square    = vp->sin_square;
        mue_square    = vp->mue_square;
        sigma_square  = vp->sigma_square;
        pos           = mobivmem[vp->vid];

        ofs.write( (char*)&p            , sizeof(p) );
        ofs.write( (char*)&radius       , sizeof(radius) );
        ofs.write( (char*)&normal       , sizeof(normal) );
        ofs.write( (char*)&sin_square   , sizeof(sin_square) );
        ofs.write( (char*)&mue_square   , sizeof(mue_square) );
        ofs.write( (char*)&sigma_square , sizeof(sigma_square) );
        ofs.write( (char*)&pos          , sizeof(pos) );

        handle2index_map[vh] = i;
    }

    for (size_t fi = 0; fi < faces_.size(); fi++) 
    {
        int fhh = halfedge_handle_fh(fi);
        int fvh = to_vertex_handle(fhh);
        fvi[0]  = handle2index_map[fvh];

        fhh     = next_halfedge_handle(fhh);
        fvh     = to_vertex_handle(fhh);
        fvi[1]  = handle2index_map[fvh];

        fhh     = next_halfedge_handle(fhh);
        fvh     = to_vertex_handle(fhh);
        fvi[2]  = handle2index_map[fvh];

        ofs.write( (char*)&fvi[0] , sizeof(fvi[0]) );
        ofs.write( (char*)&fvi[1] , sizeof(fvi[1]) );
        ofs.write( (char*)&fvi[2] , sizeof(fvi[2]) );
    }

    // write progressive detail (vertex hierarchy)

    for (size_t i = 0; i < n_details_; ++i)
    {
        VsData32 vs;

        PMInfo  pminfo = *pmiter_;

        VPIter pvp = vertex_pminfo_[pminfo.v1].vparams_ptr; 
        VPIter lvp = pvp->lchild;
        VPIter rvp = pvp->rchild;

        vs.p                = points_[pminfo.v0];

        vs.node_index       = pvp->vid;
        vs.fund_lcut_index  = pvp->fund_lcut;
        vs.fund_rcut_index  = pvp->fund_rcut;

        vs.l_radius         = lvp->radius;
        vs.l_normal         = lvp->normal;
        vs.l_sin_square     = lvp->sin_square;
        vs.l_mue_square     = lvp->mue_square;
        vs.l_sigma_square   = lvp->sigma_square;

        vs.r_radius         = rvp->radius;
        vs.r_normal         = rvp->normal;
        vs.r_sin_square     = rvp->sin_square;
        vs.r_mue_square     = rvp->mue_square;
        vs.r_sigma_square   = rvp->sigma_square;

        vs.l_pos            = mobivmem[lvp->vid];
        vs.r_pos            = mobivmem[rvp->vid];

        ofs.write( (char*)&vs , sizeof(vs) );

        refinePM(i);
    }
    ofs.close();

    std::cout << "save view-dependent progressive mesh" << std::endl;
}


void refinePM(unsigned int _n)
{
    while (n_current_res_ < _n && pmiter_ != pminfos_.end())
    {
        vertex_split_fast(pmiter_->v0,pmiter_->v1,pmiter_->vl,pmiter_->vr);

        VPIter  pvp = vertex_pminfo_[pmiter_->v1].vparams_ptr;
        vertex_pminfo_[pmiter_->v0].vparams_ptr = pvp->lchild;
        vertex_pminfo_[pmiter_->v1].vparams_ptr = pvp->rchild;

        ++pmiter_; ++n_current_res_;
    }

    // Note! We move these resize routines here to speed up a single vertex-split; (0.1.0)
    edge_status_.resize(n_edges());
    edge_pminfo_.resize(n_edges());
    face_normal_.resize(n_faces());
    face_status_.resize(n_faces());
}

void coarsePM(unsigned int _n)
{
    while (n_current_res_ > _n && pmiter_ != pminfos_.begin()) 
    {
        --pmiter_;

        int hh = find_halfedge(pmiter_->v0, pmiter_->v1); assert(to_vertex_handle(hh) == pmiter_->v1);
        collapse_fast(hh);

        VPIter  rvp = vertex_pminfo_[pmiter_->v1].vparams_ptr;
        vertex_pminfo_[pmiter_->v1].vparams_ptr = rvp->parent;

        --n_current_res_;
    }
}

inline VPIter   half_vpit (int _heh)                         { return edge_pminfo_[_heh>>1].vparams_ptr[_heh&1]                         ; }
inline void set_half_vpit (int _heh, VPIter   _vparams_ptr)  {        edge_pminfo_[_heh>>1].vparams_ptr[_heh&1]      = _vparams_ptr     ; }

void lovePM()
{
    int                             vh;
    int                             h, o, hn, op, hpo, on, ono;

    VPIterContainer                 leaf_list;


    OpenMesh::Utils::Timer tana;
    tana.start();

    refinePM(n_max_res_);

    update_normals(draw_faces_);

    std::cerr << n_vertices() << " vertices, "
              << n_edges()    << " edge, "
              << n_faces()    << " faces, "
              << n_details_   << " detail vertices\n";

    std::cout << "Init view-dependent PM analysis" << std::endl;

    // initialize
    for (size_t i = 0; i < n_halfedges(); i++)
    {
        vh = to_vertex_handle(i);

        set_half_vpit (i, vertex_pminfo_[vh].vparams_ptr);
    }

    for (size_t i = 0; i < n_vertices(); i++)
    {
        vertex_pminfo_[i].vparams_ptr->normal = vertex_normal_[i];
    }

    std::cout << "Start view-dependent PM analysis" << std::endl;

    // locate fundamental cut vertices in each edge collapse
    OpenMesh::Utils::Timer t;

    for (size_t i = n_max_res_; i>0; --i)
    {
        t.start();

        PMInfo  pminfo = pminfos_[i-1];

        if (verbose) std::cout << "Analyzing " << i << "-th detail vertex" << std::endl;

        // maintain leaf node pointers & locate fundamental cut vertices
        h   = find_halfedge(pminfo.v0, pminfo.v1);
        o   = opposite_halfedge_handle(h);
        hn  = next_halfedge_handle(h);
        hpo = opposite_halfedge_handle(prev_halfedge_handle(h));
        op  = prev_halfedge_handle(o);
        on  = next_halfedge_handle(o);
        ono = opposite_halfedge_handle(on);

        VPIter rvp = vertex_pminfo_[pminfo.v1].vparams_ptr;
        VPIter pvp = rvp->parent;

        if (pminfo.vl != -1)
        {
            VPIter fund_lcut_vpit = half_vpit(hn );
            VPIter left_leaf_vpit = half_vpit(hpo);

            set_half_vpit(hn,left_leaf_vpit);

            pvp->fund_lcut = fund_lcut_vpit->vid;
        }

        if (pminfo.vr != -1)
        {
            VPIter fund_rcut_vpit = half_vpit(on );
            VPIter right_leaf_vpit = half_vpit(ono);

            set_half_vpit(op,right_leaf_vpit);

            pvp->fund_rcut = fund_rcut_vpit->vid;
        }

        coarsePM(i-1);

        leaf_list.clear();

        get_leaf_node_handles       (pvp, leaf_list);
        compute_bounding_box        (pvp, leaf_list);
        compute_cone_of_normals     (pvp, leaf_list);
        compute_screen_space_error  (pvp, leaf_list);

        t.stop();

        if (verbose)
        {
            std::cout << "  radius of bounding sphere: "                << pvp->radius << std::endl;
            std::cout << "  direction of cone of normals: "             << pvp->normal << std::endl;
            std::cout << "  sin(semi-angle of cone of normals) ^2: "    << pvp->sin_square << std::endl;
            std::cout << "  (mue^2, sigma^2) : ("                       << pvp->mue_square   << ", " 
                                                                        << pvp->sigma_square << ")"  << std::endl;
            std::cout << "- " << t.as_string() << std::endl;
        }

    } // end for all collapses

    tana.stop();

    std::cout << "Analyzing step completed in " << tana.as_string() << std::endl;
}



void get_leaf_node_handles (VPIter vp, VPIterContainer &leaf_list)
{
    if (vp->lchild == vplist.end() && vp->rchild == vplist.end())
    {
        leaf_list.push_back(vp);
    }
    else
    {
        get_leaf_node_handles(vp->lchild, leaf_list);
        get_leaf_node_handles(vp->rchild, leaf_list);
    }
}

void compute_bounding_box (VPIter vp, VPIterContainer &leaf_list)
{
    float                       max_distance;
    Vec3f                       p, lp;

    max_distance = 0.0f; int vh = vp->vh; p = points_[vh];

    for (size_t i = 0; i < leaf_list.size(); i++ )
    {
        lp              = points_[leaf_list[i]->vh];
        max_distance    = std::max(max_distance, (p - lp).length());
    }
    vp->radius = max_distance;  
}

void compute_cone_of_normals (VPIter vp, VPIterContainer &leaf_list)
{
    float                           max_angle, angle;
    Vec3f                           n, ln;
    int                             vh = vp->vh;

    // Note! We calculate the vertex normal here; (0.1.0)
    n = Vec3f(0.0f,0.0f,0.0f); unsigned int nsize = 0;
        int shh = vertices_[vh].halfedge_handle_;
        int  hh = shh;
        do {
            const Edge& nedge = edges_[hh >> 1];
            int fh = nedge.halfedges_[hh & 1].face_handle_; 

            const Vec3f&    p0 = points_[draw_faces_[3*fh  ]]; 
            const Vec3f&    p1 = points_[draw_faces_[3*fh+1]];
            const Vec3f&    p2 = points_[draw_faces_[3*fh+2]]; 
            Vec3f         p1p0 = p0 - p1;
            Vec3f         p1p2 = p2 - p1;
            Vec3f          tfn = p1p2 % p1p0; float norm = tfn.length();
            Vec3f           fn = ( (norm != float(0)) ? ((tfn *= (float(1)/norm)),tfn) : Vec3f(0,0,0) );

            n += fn; nsize++;

            hh = nedge.halfedges_[(hh&1) ^ 1].next_halfedge_handle_;
        } while (hh != shh);
    assert(nsize!=0); n = n/nsize;

    max_angle = 0.0f;

    for (size_t i = 0; i < leaf_list.size(); i++ )
    {
        ln        = leaf_list[i]->normal;
        angle     = acosf( dot(n,ln) );
        max_angle = std::max(max_angle, angle );
    }

    vertex_normal_[vh]  = n;
    vp->normal          = n;

    max_angle           = std::min(max_angle, float(M_PI_2));
    vp->sin_square      = sinf(max_angle)*sinf(max_angle);
}

void compute_screen_space_error (VPIter vp, VPIterContainer &leaf_list)
{
    std::vector<Vec3f>      residuals;
    Vec3f                   residual, res;
    Vec3f                   lp, tri[3];
    float                   min_distance;
    float                   s, t;

    int            vh = vp->vh;
    const Vec3f&    p = points_[vh];

    std::vector<Vec3f> vfinfos;
    int shh = vertices_[vh].halfedge_handle_;
    int  hh = shh;
    do {
        const Edge& nedge = edges_[hh >> 1];
        int fh = nedge.halfedges_[hh & 1].face_handle_; 

        tri[0] = points_[draw_faces_[3*fh  ]]; 
        tri[1] = points_[draw_faces_[3*fh+1]];
        tri[2] = points_[draw_faces_[3*fh+2]]; 

        Vec3f    B = tri[0];            vfinfos.push_back(B );  // Tri.Origin();
        Vec3f   E0 = tri[1] - tri[0];   vfinfos.push_back(E0);  // rkTri.Edge0()
        Vec3f   E1 = tri[2] - tri[0];   vfinfos.push_back(E1);  // rkTri.Edge1()

        hh = nedge.halfedges_[(hh&1) ^ 1].next_halfedge_handle_;
    } while (hh != shh);

    for (size_t i = 0; i < leaf_list.size(); i++)
    {
        lp = points_[leaf_list[i]->vh];

        // compute residual of a leaf-vertex from the current mesh_
        residual = lp - p; min_distance = residual.length();

        for (size_t j = 0; j < vfinfos.size(); j+=3)
        {
            res = point2triangle_residual(lp, vfinfos[j], vfinfos[j+1], vfinfos[j+2], s, t);

            if (res.length() < min_distance) {residual = res;min_distance = res.length();}
        }

        residuals.push_back(residual);
    }
    compute_mue_sigma(vp, residuals);
}

void compute_mue_sigma(VPIter vp, ResidualContainer &residuals)
{
    Vec3f       vn;
    float       max_inner, max_cross;

    ResidualContainer::iterator  r_it, r_end(residuals.end());

    max_inner = max_cross = 0.0f;
    vn = vertex_normal_[vp->vh];

    for (r_it = residuals.begin(); r_it != r_end; ++r_it)
    {
        float inner = fabsf(dot(*r_it, vn));
        float cross = OpenMesh::cross(*r_it, vn).length();

        max_inner = std::max(max_inner, inner);
        max_cross = std::max(max_cross, cross);
    }
    if (max_cross < 1.0e-7)
    {
        vp->mue_square   = max_cross * max_cross;
        vp->sigma_square = max_inner * max_inner;
    }
    else 
    {
        float  ratio = std::max(1.0f, max_inner/max_cross);
        float  whole_degree = acosf(1.0f/ratio);
        float  mue, max_mue;
        float  degree;
        float  res_length;
        Vec3f  res;

        max_mue = 0.0f;
        for (r_it = residuals.begin(); r_it != r_end; ++r_it)
        {
            res = *r_it;
            res_length = res.length();

            // TODO: take care when res.length() is too small
            degree = acosf(dot(vn,res) / res_length);

            if (degree < 0.0f         ) degree =            -degree;
            if (degree > float(M_PI_2)) degree = float(M_PI)-degree;

            mue     = (degree < whole_degree) ? cosf(whole_degree - degree) * res_length : res_length;

            max_mue = std::max(max_mue, mue); 
        }

        vp->mue_square   = max_mue * max_mue;
        vp->sigma_square = (ratio*max_mue) * (ratio*max_mue);
    }
}

Vec3f point2triangle_residual(const Vec3f &p, const Vec3f& B, const Vec3f& E0, const Vec3f& E1, float &s, float &t)
{ 
    // OpenMesh::Vec3f B  = tri[0];                 // Tri.Origin();
    // OpenMesh::Vec3f E0 = tri[1] - tri[0];        // rkTri.Edge0()
    // OpenMesh::Vec3f E1 = tri[2] - tri[0];        // rkTri.Edge1()
    // OpenMesh::Vec3f D  = tri[0] - p;             // kDiff

    Vec3f   D = B - p;                  // kDiff

    float	a = dot(E0, E0);                    // fA00
    float	b = dot(E0, E1);                    // fA01
    float	c = dot(E1, E1);                    // fA11
    float	d = dot(E0, D );                    // fB0
    float	e = dot(E1, D );                    // fB1
    float	f = dot(D , D );                    // fC
    float det = fabsf(a*c - b*b);
    s = b*e-c*d;
    t = b*d-a*e;

    OpenMesh::Vec3f  residual;

    float distance2;

    if ( s + t <= det )
    {
        if ( s < 0.0f )
        {
            if ( t < 0.0f )  // region 4
            {
                if ( d < 0.0f )
                {
                    t = 0.0f;
                    if ( -d >= a ) {s = 1.0f; distance2 = a+2.0f*d+f;}
                    else           {s = -d/a; distance2 = d*s+f     ;}
                } 
                else
                {
                    s = 0.0f;
                    if      (  e >= 0.0f ) {t = 0.0f; distance2 = f         ;}
                    else if ( -e >= c    ) {t = 1.0f; distance2 = c+2.0f*e+f;}
                    else                   {t = -e/c; distance2 = e*t+f     ;}
                }
            }
            else  // region 3
            {
                s = 0.0f;
                if      (  e >= 0.0f ) {t = 0.0f; distance2 = f         ;}
                else if ( -e >= c    ) {t = 1.0f; distance2 = c+2.0f*e+f;}
                else                   {t = -e/c; distance2 = e*t+f     ;}
            }
        }
        else if ( t < 0.0f )  // region 5
        {
            t = 0.0f;
            if ( d >= 0.0f ){s = 0.0f;distance2 = f;}
            else if ( -d >= a ){s = 1.0f;distance2 = a+2.0f*d+f; }
            else {s = -d/a;distance2 = d*s+f;}
        }
        else  // region 0
        {
            // minimum at interior point
            float inv_det = 1.0f/det;s *= inv_det;t *= inv_det;
            distance2 = s*(a*s+b*t+2.0f*d) + t*(b*s+c*t+2.0f*e)+f;
        }
    }
    else
    {
        float tmp0, tmp1, numer, denom;

        if ( s < 0.0f )  // region 2
        {
            tmp0 = b + d;
            tmp1 = c + e;
            if ( tmp1 > tmp0 )
            {
                numer = tmp1 - tmp0;
                denom = a-2.0f*b+c;
                if ( numer >= denom ) {s = 1.0f       ; t = 0.0f    ; distance2 = a+2.0f*d+f                              ;}
                else                  {s = numer/denom; t = 1.0f - s; distance2 = s*(a*s+b*t+2.0f*d) +t*(b*s+c*t+2.0f*e)+f;}
            }
            else
            {
                s = 0.0f;
                if   ( tmp1 <= 0.0f ) {t = 1.0f; distance2 = c+2.0f*e+f;}
                else if ( e >= 0.0f ) {t = 0.0f; distance2 = f         ;}
                else                  {t = -e/c; distance2 = e*t+f     ;}
            }
        }
        else if ( t < 0.0f )  // region 6
        {
            tmp0 = b + e;
            tmp1 = a + d;
            if ( tmp1 > tmp0 )
            {
                numer = tmp1 - tmp0;
                denom = a-2.0f*b+c;
                if ( numer >= denom ) {t = 1.0f       ; s = 0.0f    ; distance2 = c+2.0f*e+f                              ;}
                else                  {t = numer/denom; s = 1.0f - t; distance2 = s*(a*s+b*t+2.0f*d)+ t*(b*s+c*t+2.0f*e)+f;}
            }
            else
            {
                t = 0.0f;
                if      ( tmp1 <= 0.0f ) {s = 1.0f; distance2 = a+2.0f*d+f;}
                else if (    d >= 0.0f ) {s = 0.0f; distance2 = f         ;}
                else                     {s = -d/a; distance2 = d*s+f     ;}
            }
        }
        else  // region 1
        {
            numer = c + e - b - d;
            if     ( numer <= 0.0f  ) {s = 0.0f       ; t = 1.0f    ; distance2 = c+2.0f*e+f                               ;}
            else {
                denom = a-2.0f*b+c;
                if ( numer >= denom ) {s = 1.0f       ; t = 0.0f    ; distance2 = a+2.0f*d+f                               ;}
                else                  {s = numer/denom; t = 1.0f - s; distance2 = s*(a*s+b*t+2.0f*d) + t*(b*s+c*t+2.0f*e)+f;}
            }
        }
    }

    residual = p - (B + s*E0 + t*E1);

    return	residual;
}

/* Vec3f point2triangle_residual(const Vec3f &p, const Vec3f tri[3], float &s, float &t)
{ 
    OpenMesh::Vec3f B = tri[0];                 // Tri.Origin();
    OpenMesh::Vec3f E0 = tri[1] - tri[0];       // rkTri.Edge0()
    OpenMesh::Vec3f E1 = tri[2] - tri[0];       // rkTri.Edge1()
    OpenMesh::Vec3f D = tri[0] - p;             // kDiff
    float	a = dot(E0, E0);                    // fA00
    float	b = dot(E0, E1);                    // fA01
    float	c = dot(E1, E1);                    // fA11
    float	d = dot(E0, D);                     // fB0
    float	e = dot(E1, D);                     // fB1
    float	f = dot(D, D);                      // fC
    float det = fabsf(a*c - b*b);
    s = b*e-c*d;
    t = b*d-a*e;

    OpenMesh::Vec3f     residual;

    float distance2;

    if ( s + t <= det )
    {
        if ( s < 0.0f )
        {
            if ( t < 0.0f )  // region 4
            {
                if ( d < 0.0f )
                {
                    t = 0.0f;
                    if ( -d >= a ) {s = 1.0f; distance2 = a+2.0f*d+f;}
                    else           {s = -d/a; distance2 = d*s+f     ;}
                } 
                else
                {
                    s = 0.0f;
                    if      (  e >= 0.0f ) {t = 0.0f; distance2 = f         ;}
                    else if ( -e >= c    ) {t = 1.0f; distance2 = c+2.0f*e+f;}
                    else                   {t = -e/c; distance2 = e*t+f     ;}
                }
            }
            else  // region 3
            {
                s = 0.0f;
                if      (  e >= 0.0f ) {t = 0.0f; distance2 = f         ;}
                else if ( -e >= c    ) {t = 1.0f; distance2 = c+2.0f*e+f;}
                else                   {t = -e/c; distance2 = e*t+f     ;}
            }
        }
        else if ( t < 0.0f )  // region 5
        {
            t = 0.0f;
            if ( d >= 0.0f ){s = 0.0f;distance2 = f;}
            else if ( -d >= a ){s = 1.0f;distance2 = a+2.0f*d+f; }
            else {s = -d/a;distance2 = d*s+f;}
        }
        else  // region 0
        {
            // minimum at interior point
            float inv_det = 1.0f/det;s *= inv_det;t *= inv_det;
            distance2 = s*(a*s+b*t+2.0f*d) + t*(b*s+c*t+2.0f*e)+f;
        }
    }
    else
    {
        float tmp0, tmp1, numer, denom;

        if ( s < 0.0f )  // region 2
        {
            tmp0 = b + d;
            tmp1 = c + e;
            if ( tmp1 > tmp0 )
            {
                numer = tmp1 - tmp0;
                denom = a-2.0f*b+c;
                if ( numer >= denom ) {s = 1.0f       ; t = 0.0f    ; distance2 = a+2.0f*d+f                              ;}
                else                  {s = numer/denom; t = 1.0f - s; distance2 = s*(a*s+b*t+2.0f*d) +t*(b*s+c*t+2.0f*e)+f;}
            }
            else
            {
                s = 0.0f;
                if   ( tmp1 <= 0.0f ) {t = 1.0f; distance2 = c+2.0f*e+f;}
                else if ( e >= 0.0f ) {t = 0.0f; distance2 = f         ;}
                else                  {t = -e/c; distance2 = e*t+f     ;}
            }
        }
        else if ( t < 0.0f )  // region 6
        {
            tmp0 = b + e;
            tmp1 = a + d;
            if ( tmp1 > tmp0 )
            {
                numer = tmp1 - tmp0;
                denom = a-2.0f*b+c;
                if ( numer >= denom ) {t = 1.0f       ; s = 0.0f    ; distance2 = c+2.0f*e+f                              ;}
                else                  {t = numer/denom; s = 1.0f - t; distance2 = s*(a*s+b*t+2.0f*d)+ t*(b*s+c*t+2.0f*e)+f;}
            }
            else
            {
                t = 0.0f;
                if      ( tmp1 <= 0.0f ) {s = 1.0f; distance2 = a+2.0f*d+f;}
                else if (    d >= 0.0f ) {s = 0.0f; distance2 = f         ;}
                else                     {s = -d/a; distance2 = d*s+f     ;}
            }
        }
        else  // region 1
        {
            numer = c + e - b - d;
            if     ( numer <= 0.0f  ) {s = 0.0f       ; t = 1.0f    ; distance2 = c+2.0f*e+f                               ;}
            else {
                denom = a-2.0f*b+c;
                if ( numer >= denom ) {s = 1.0f       ; t = 0.0f    ; distance2 = a+2.0f*d+f                               ;}
                else                  {s = numer/denom; t = 1.0f - s; distance2 = s*(a*s+b*t+2.0f*d) + t*(b*s+c*t+2.0f*e)+f;}
            }
        }
    }

    residual = p - (B + s*E0 + t*E1);

    return	residual;
} */

// ============================================================================
