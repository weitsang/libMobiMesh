#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <sstream>
#include <cmath>
#include <iomanip>
#include <cstdlib>
#include <fcntl.h>
#include <unistd.h>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <errno.h>

#include <vector>
#include <list>
#include <bitset>
#include <set>
#include <map>
#include <queue>
#include <stack>
#include <string>
#include <cstring>
#include <algorithm>

// Note! We define GL_GLEXT_PROTOTYPES to get function prototypes from glext.h; (0.6.0)
//#define GL_GLEXT_PROTOTYPES

//#include <GL/glut.h>
//#include <GL/glext.h>

// Note! We remove most OpenMesh features except for VectorT; (0.4.2)
#include "inc/VectorT.hh"

using namespace OpenMesh;

class Endian {
public:
    enum Type {
        LSB = 1, // Little endian (Intel family and clones)
        MSB      // Big endian (Motorola's 68x family, DEC Alpha, MIPS)
    };
    static Type local() { int one_ = 1; return *((unsigned char*)&one_) ? Endian::LSB : Endian::MSB; }
};
/*
class VHierarchyNodeIndex {

private: unsigned int value_;
public:  static const VHierarchyNodeIndex InvalidIndex;
public:
    VHierarchyNodeIndex() { value_ = 0; }
    VHierarchyNodeIndex(unsigned int _value) { value_ = _value; }
    VHierarchyNodeIndex(const VHierarchyNodeIndex &_other) { value_ = _other.value_; }
    VHierarchyNodeIndex(unsigned int _tree_id, unsigned int _node_id, unsigned short _tree_id_bits) {
        assert(_tree_id < ((unsigned int) 0x00000001 << _tree_id_bits));
        assert(_node_id < ((unsigned int) 0x00000001 << (32 - _tree_id_bits)));
        value_ = (_tree_id << (32 - _tree_id_bits)) | _node_id; 
    }
    bool is_valid(unsigned short _tree_id_bits)        const { return  node_id(_tree_id_bits) != 0 ? true : false;  }
    unsigned int tree_id(unsigned short _tree_id_bits) const { return  value_ >> (32 - _tree_id_bits); }
    unsigned int node_id(unsigned short _tree_id_bits) const { return  value_ & ((unsigned int) 0xFFFFFFFF >> _tree_id_bits); }
    bool operator< (const VHierarchyNodeIndex &other)  const { return  (value_ < other.value_) ? true : false; }
    unsigned int value() const { return  value_; }
};
const VHierarchyNodeIndex VHierarchyNodeIndex::InvalidIndex = VHierarchyNodeIndex();
*/
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

typedef struct VsData64 {
    Vec3f                       p;
    unsigned int                padding;
    unsigned long long          node_index;
    unsigned long long          fund_lcut_index; //  left cut neighbor
    unsigned long long          fund_rcut_index; // right cut neighbor
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
} VsData64;

// Note! We load view-dependent parameters from secondary storage when necessary;
typedef struct ZsData32 {
    Vec3f                       p;
    unsigned int                fund_lcut_index;
    unsigned int                fund_rcut_index;
    Vec3f                       l_normal;
    Vec3f                       r_normal;
    unsigned int                l_pos;
    unsigned int                r_pos;
} ZsData32;

typedef struct ZsData64 {
    Vec3f                       p;
    unsigned int                padding;
    unsigned long long          fund_lcut_index;
    unsigned long long          fund_rcut_index;
    Vec3f                       l_normal;
    Vec3f                       r_normal;
    unsigned int                l_pos;
    unsigned int                r_pos;
} ZsData64;

typedef struct ZipZsData32 {
    Vec3f                       p;
    unsigned int                fund_lcut_index;
    unsigned int                fund_rcut_index;
    int                         c_pos;
} ZipZsData32;

typedef struct ZipZsData64 {
    Vec3f                       p;
    unsigned int                fund_ccut_tree;
    unsigned int                fund_lcut_node;
    unsigned int                fund_rcut_node;
    int                         c_pos;
} ZipZsData64;

typedef struct ZipTVData32 {
    Vec3f                       p;
    Vec3f                       l_normal;
    Vec3f                       r_normal;
    unsigned int                fund_lcut_index;
    unsigned int                fund_rcut_index;
    int                         c_pos;
} ZipTVData32;

typedef struct ZipTVData64 {
    Vec3f                       p;
    Vec3f                       l_normal;
    Vec3f                       r_normal;
    unsigned int                fund_ccut_tree;
    unsigned int                fund_lcut_node;
    unsigned int                fund_rcut_node;
    int                         c_pos;
} ZipTVData64;

// Note! We use 'unsigned int' (32-bit) or 'unsigned long (long)' (64-bit) as VertexID; (0.5.4)
typedef unsigned int            VertexID32;
typedef unsigned long long      VertexID64;

void VPM322ZPM32 (const char * ifile, const char * ofile)
{
    unsigned int        n_base_vertices_;
    unsigned int        n_base_faces_;
    unsigned int        n_details_;

    char                fileformat[16];
    unsigned int        fvi[3];
    unsigned int        pos;
    Vec3f               p, normal;
    float               radius, sin_square, mue_square, sigma_square;

    std::ifstream ifs(ifile,    std::ios::binary);
    std::ofstream ofs(ofile,    std::ios::binary); 

    if (!ifs) { std::cerr << "read error\n"; exit(1); }
    ifs.read(fileformat, 10); fileformat[10] = '\0';
    if (std::string(fileformat) != std::string("VDProgMesh"))
    { std::cerr << "Wrong file format.\n"; ifs.close(); exit(1); }

    ofs.write(fileformat, 10); // ofs << "VDProgMesh";

    ifs. read( (char*)&n_base_vertices_, sizeof(n_base_vertices_) );
    ifs. read( (char*)&n_base_faces_   , sizeof(n_base_faces_   ) );
    ifs. read( (char*)&n_details_      , sizeof(n_details_      ) );

    ofs.write( (char*)&n_base_vertices_, sizeof(n_base_vertices_) );
    ofs.write( (char*)&n_base_faces_   , sizeof(n_base_faces_   ) );
    ofs.write( (char*)&n_details_      , sizeof(n_details_      ) );

    for (size_t i=0; i<n_base_vertices_; ++i)
    {
        ifs. read( (char*)&p            , sizeof(p) );
        ifs. read( (char*)&radius       , sizeof(radius) );
        ifs. read( (char*)&normal       , sizeof(normal) );
        ifs. read( (char*)&sin_square   , sizeof(sin_square) );
        ifs. read( (char*)&mue_square   , sizeof(mue_square) );
        ifs. read( (char*)&sigma_square , sizeof(sigma_square) );
        ifs. read( (char*)&pos          , sizeof(pos) );

        ofs.write( (char*)&p            , sizeof(p) );
        ofs.write( (char*)&normal       , sizeof(normal) );
        ofs.write( (char*)&pos          , sizeof(pos) );
    }

    for (size_t i=0; i<n_base_faces_; ++i)
    {
        ifs.read ( (char*)&fvi[0] , sizeof(fvi[0]) );
        ifs.read ( (char*)&fvi[1] , sizeof(fvi[1]) );
        ifs.read ( (char*)&fvi[2] , sizeof(fvi[2]) );

        ofs.write( (char*)&fvi[0] , sizeof(fvi[0]) );
        ofs.write( (char*)&fvi[1] , sizeof(fvi[1]) );
        ofs.write( (char*)&fvi[2] , sizeof(fvi[2]) );
    }
    for (size_t i=0; i<n_details_; i++)
    {
        VsData32 vs; ifs.read( (char*)&vs , sizeof(vs) );
        ZsData32 zs; 
        zs.              p = vs.p;
        zs.fund_lcut_index = vs.fund_lcut_index; //  left cut neighbor
        zs.fund_rcut_index = vs.fund_rcut_index; // right cut neighbor
        zs.       l_normal = vs.l_normal;
        zs.       r_normal = vs.r_normal;
        zs.          l_pos = vs.l_pos;
        zs.          r_pos = vs.r_pos;
        ofs.write( (char*)&zs , sizeof(zs) );
    }
    ifs.close(); ofs.close();

    std::cerr << n_base_vertices_ << " vertices, "
              << n_base_faces_    << " faces, "
              << n_details_       << " detail vertices\n";
}

void VPM642ZPM64 (const char * ifile, const char * ofile) 
{
    unsigned long long          n_base_vertices_;
    unsigned long long          n_base_faces_;
    unsigned long long          n_details_;

    char                        fileformat[16];
    unsigned long long          fvi[3];
    unsigned int                pos;
    Vec3f                       p, normal;
    float                       radius, sin_square, mue_square, sigma_square;

    std::ifstream ifs(ifile,    std::ios::binary);
    std::ofstream ofs(ofile,    std::ios::binary); 

    if (!ifs) { std::cerr << "read error\n"; exit(1); }
    ifs.read(fileformat, 10); fileformat[10] = '\0';
    if (std::string(fileformat) != std::string("VDProgMesh"))
    { std::cerr << "Wrong file format.\n"; ifs.close(); exit(1); }

    ofs.write(fileformat, 10); // ofs << "VDProgMesh";

    ifs. read( (char*)&n_base_vertices_, sizeof(n_base_vertices_) );
    ifs. read( (char*)&n_base_faces_   , sizeof(n_base_faces_   ) );
    ifs. read( (char*)&n_details_      , sizeof(n_details_      ) );

    ofs.write( (char*)&n_base_vertices_, sizeof(n_base_vertices_) );
    ofs.write( (char*)&n_base_faces_   , sizeof(n_base_faces_   ) );
    ofs.write( (char*)&n_details_      , sizeof(n_details_      ) );

    for (size_t i=0; i<n_base_vertices_; ++i)
    {
        ifs. read( (char*)&p            , sizeof(p) );
        ifs. read( (char*)&radius       , sizeof(radius) );
        ifs. read( (char*)&normal       , sizeof(normal) );
        ifs. read( (char*)&sin_square   , sizeof(sin_square) );
        ifs. read( (char*)&mue_square   , sizeof(mue_square) );
        ifs. read( (char*)&sigma_square , sizeof(sigma_square) );
        ifs. read( (char*)&pos          , sizeof(pos) );

        ofs.write( (char*)&p            , sizeof(p) );
        ofs.write( (char*)&normal       , sizeof(normal) );
        ofs.write( (char*)&pos          , sizeof(pos) );
    }

    for (size_t i=0; i<n_base_faces_; ++i)
    {
        ifs.read ( (char*)&fvi[0] , sizeof(fvi[0]) );
        ifs.read ( (char*)&fvi[1] , sizeof(fvi[1]) );
        ifs.read ( (char*)&fvi[2] , sizeof(fvi[2]) );

        ofs.write( (char*)&fvi[0] , sizeof(fvi[0]) );
        ofs.write( (char*)&fvi[1] , sizeof(fvi[1]) );
        ofs.write( (char*)&fvi[2] , sizeof(fvi[2]) );
    }
    for (size_t i=0; i<n_details_; i++)
    {
        VsData64 vs; ifs.read( (char*)&vs , sizeof(vs) );
        ZsData64 zs; 
        zs.              p = vs.p;
        zs.        padding = vs.padding;
        zs.fund_lcut_index = vs.fund_lcut_index; //  left cut neighbor
        zs.fund_rcut_index = vs.fund_rcut_index; // right cut neighbor
        zs.       l_normal = vs.l_normal;
        zs.       r_normal = vs.r_normal;
        zs.          l_pos = vs.l_pos;
        zs.          r_pos = vs.r_pos;
        ofs.write( (char*)&zs , sizeof(zs) );
    }
    ifs.close(); ofs.close();

    std::cerr << n_base_vertices_ << " vertices, "
              << n_base_faces_    << " faces, "
              << n_details_       << " detail vertices\n";
}

void ZipVPM322ZPM32 (const char * ifile, const char * ofile)
{
    unsigned int        n_base_vertices_;
    unsigned int        n_base_faces_;
    unsigned int        n_details_;

    char                fileformat[16];
    unsigned int        fvi[3];
    unsigned int        pos;
    Vec3f               p, normal;
    float               radius, sin_square, mue_square, sigma_square;

    std::ifstream ifs(ifile,    std::ios::binary);
    std::ofstream ofs(ofile,    std::ios::binary); 

    if (!ifs) { std::cerr << "read error\n"; exit(1); }
    ifs.read(fileformat, 10); fileformat[10] = '\0';
    if (std::string(fileformat) != std::string("VDProgMesh"))
    { std::cerr << "Wrong file format.\n"; ifs.close(); exit(1); }

    ofs.write(fileformat, 10); // ofs << "VDProgMesh";

    ifs. read( (char*)&n_base_vertices_, sizeof(n_base_vertices_) );
    ifs. read( (char*)&n_base_faces_   , sizeof(n_base_faces_   ) );
    ifs. read( (char*)&n_details_      , sizeof(n_details_      ) );

    ofs.write( (char*)&n_base_vertices_, sizeof(n_base_vertices_) );
    ofs.write( (char*)&n_base_faces_   , sizeof(n_base_faces_   ) );
    ofs.write( (char*)&n_details_      , sizeof(n_details_      ) );

    size_t s_params_, s_base_vertices_, s_base_faces_, s_base_, s_details_;

    s_params_        = sizeof("VDProgMesh")-1+sizeof(unsigned int)+sizeof(unsigned int)+sizeof(unsigned int);
    s_base_vertices_ = sizeof(Vec3f)+sizeof(float)+sizeof(Vec3f)+sizeof(float)+sizeof(float)+sizeof(float)+sizeof(unsigned int);
    s_base_faces_    = sizeof(unsigned int)+sizeof(unsigned int)+sizeof(unsigned int);

    s_base_          = s_params_+n_base_vertices_*s_base_vertices_+n_base_faces_*s_base_faces_;
    s_details_       = sizeof(Vec3f)+sizeof(unsigned int)+sizeof(unsigned int)+sizeof(unsigned int) 
                      +sizeof(float)+sizeof(Vec3f)+sizeof(float)+sizeof(float)+sizeof(float)
                      +sizeof(float)+sizeof(Vec3f)+sizeof(float)+sizeof(float)+sizeof(float)
                      +sizeof(unsigned int)+sizeof(unsigned int);

    assert (s_details_ == sizeof(VsData32));

    // Note! We use mmap in Linux to map files and manage caches (out-of-core); (0.1.0)
    int fd; char *data; struct stat sbuf;

    if ((fd = open(ifile, O_RDONLY)) == -1)                                                     {perror("openfile");exit(1);}
    if (stat(ifile, &sbuf) == -1)                                                               {perror("filesize");exit(1);}
    if ((data = (char *) mmap(NULL, sbuf.st_size, PROT_READ, MAP_SHARED, fd, 0)) == MAP_FAILED) {perror("mmapfile");exit(1);}

    std::vector <unsigned int>  pos_vec;
    std::vector <unsigned int>  pos_map (n_base_vertices_+n_details_+1);

    for (size_t i=0; i<n_base_vertices_; ++i)
    {
        ifs. read( (char*)&p            , sizeof(p) );
        ifs. read( (char*)&radius       , sizeof(radius) );
        ifs. read( (char*)&normal       , sizeof(normal) );
        ifs. read( (char*)&sin_square   , sizeof(sin_square) );
        ifs. read( (char*)&mue_square   , sizeof(mue_square) );
        ifs. read( (char*)&sigma_square , sizeof(sigma_square) );
        ifs. read( (char*)&pos          , sizeof(pos) );

        unsigned int zpos = 0; if (pos != 0) { pos_vec.push_back(pos); pos_map[pos] = pos_vec.size(); zpos = pos_map[pos]; }

        ofs.write( (char*)&p            , sizeof(p) );
        ofs.write( (char*)&zpos         , sizeof(zpos) );
    }

    for (size_t i=0; i<n_base_faces_; ++i)
    {
        ifs.read ( (char*)&fvi[0] , sizeof(fvi[0]) );
        ifs.read ( (char*)&fvi[1] , sizeof(fvi[1]) );
        ifs.read ( (char*)&fvi[2] , sizeof(fvi[2]) );

        ofs.write( (char*)&fvi[0] , sizeof(fvi[0]) );
        ofs.write( (char*)&fvi[1] , sizeof(fvi[1]) );
        ofs.write( (char*)&fvi[2] , sizeof(fvi[2]) );
    }

    for (size_t i=0; i<n_details_; i++)
    {
        VsData32 vs; ifs.read( (char*)&vs , sizeof(vs) ); assert (!ifs.fail()); 

        if (vs.l_pos != 0) { pos_vec.push_back(vs.l_pos); pos_map[vs.l_pos] = pos_vec.size(); }
        if (vs.r_pos != 0) { pos_vec.push_back(vs.r_pos); pos_map[vs.r_pos] = pos_vec.size(); } 
    }

    for (size_t i=0; i<pos_vec.size(); i++)
    {
        const VsData32& vs = (VsData32 &) * (data + s_base_ + (pos_vec[i]-1)*s_details_);

        int zpos = 0;
        if      (vs.l_pos != 0 && vs.r_pos != 0) {zpos =  pos_map[vs.l_pos]*2+1;}
        else if (vs.l_pos != 0 && vs.r_pos == 0) {zpos =  pos_map[vs.l_pos]*2  ;}
        else if (vs.l_pos == 0 && vs.r_pos != 0) {zpos = -pos_map[vs.r_pos]*2  ;}

        ZipZsData32 zs; 
        zs.             p = vs.p;
        zs.fund_lcut_index = vs.fund_lcut_index;
        zs.fund_rcut_index = vs.fund_rcut_index;
        zs.          c_pos = zpos;
        // zs.       l_pos = vs.l_pos != 0 ? pos_map[vs.l_pos] : 0;
        // zs.       r_pos = vs.r_pos != 0 ? pos_map[vs.r_pos] : 0;

        ofs.write( (char*)&zs , sizeof(zs) );
    }
    if (munmap(data, sbuf.st_size) == -1) {perror("Error un-mmapping");} close(fd);

    ifs.close(); ofs.close();

    std::cerr << n_base_vertices_ << " vertices, "
              << n_base_faces_    << " faces, "
              << n_details_       << " detail vertices\n";
}

void ZipVPM642ZPM64 (const char * ifile, const char * ofile)  
{
    unsigned long long          n_base_vertices_;
    unsigned long long          n_base_faces_;
    unsigned long long          n_details_;

    char                        fileformat[16];
    unsigned long long          fvi[3];
    unsigned int                pos;
    Vec3f                       p, normal;
    float                       radius, sin_square, mue_square, sigma_square;

    std::ifstream ifs(ifile,    std::ios::binary);
    std::ofstream ofs(ofile,    std::ios::binary); 

    if (!ifs) { std::cerr << "read error\n"; exit(1); }
    ifs.read(fileformat, 10); fileformat[10] = '\0';
    if (std::string(fileformat) != std::string("VDProgMesh"))
    { std::cerr << "Wrong file format.\n"; ifs.close(); exit(1); }

    ofs.write(fileformat, 10); // ofs << "VDProgMesh";

    ifs. read( (char*)&n_base_vertices_, sizeof(n_base_vertices_) );
    ifs. read( (char*)&n_base_faces_   , sizeof(n_base_faces_   ) );
    ifs. read( (char*)&n_details_      , sizeof(n_details_      ) );

    unsigned int n_base_vertices =  n_base_vertices_;
    unsigned int n_base_faces    =  n_base_faces_;
    unsigned int n_details       =  n_details_;

    ofs.write( (char*)&n_base_vertices , sizeof(n_base_vertices ) );
    ofs.write( (char*)&n_base_faces    , sizeof(n_base_faces    ) );
    ofs.write( (char*)&n_details       , sizeof(n_details       ) );

    size_t s_params_, s_base_vertices_, s_base_faces_, s_base_, s_details_;

    s_params_        = sizeof("VDProgMesh")-1+sizeof(unsigned long long)+sizeof(unsigned long long)+sizeof(unsigned long long);
    s_base_vertices_ = sizeof(Vec3f)+sizeof(float)+sizeof(Vec3f)+sizeof(float)+sizeof(float)+sizeof(float)+sizeof(unsigned int);
    s_base_faces_    = sizeof(unsigned long long)+sizeof(unsigned long long)+sizeof(unsigned long long);

    s_base_          = s_params_+n_base_vertices_*s_base_vertices_+n_base_faces_*s_base_faces_;
    s_details_       = sizeof(Vec3f)+sizeof(unsigned int)+sizeof(unsigned long long)+sizeof(unsigned long long)+sizeof(unsigned long long) 
                      +sizeof(float)+sizeof(Vec3f)+sizeof(float)+sizeof(float)+sizeof(float)
                      +sizeof(float)+sizeof(Vec3f)+sizeof(float)+sizeof(float)+sizeof(float)
                      +sizeof(unsigned int)+sizeof(unsigned int);

    assert (s_details_ == sizeof(VsData64)); 

    // Note! We use mmap in Linux to map files and manage caches (out-of-core); (0.1.0)
    int fd; char *data; struct stat sbuf;

    if ((fd = open(ifile, O_RDONLY)) == -1)                                                     {perror("openfile");exit(1);}
    if (stat(ifile, &sbuf) == -1)                                                               {perror("filesize");exit(1);}
    if ((data = (char *) mmap(NULL, sbuf.st_size, PROT_READ, MAP_SHARED, fd, 0)) == MAP_FAILED) {perror("mmapfile");exit(1);}

    unsigned char       tree_id_bits_     = 0; while(n_base_vertices_ > ((unsigned long long) 1 << tree_id_bits_)) ++tree_id_bits_;
    unsigned char       tree_id_bitshift_ = sizeof(unsigned long long)*8   - tree_id_bits_;
    unsigned long long  tree_id_mask_     =      ((unsigned long long)-1) >> tree_id_bits_;

    std::vector <unsigned int>  pos_vec;
    std::vector <unsigned int>  pos_map (n_base_vertices_+n_details_+1);

    for (size_t i=0; i<n_base_vertices_; ++i)
    {
        ifs. read( (char*)&p            , sizeof(p) );
        ifs. read( (char*)&radius       , sizeof(radius) );
        ifs. read( (char*)&normal       , sizeof(normal) );
        ifs. read( (char*)&sin_square   , sizeof(sin_square) );
        ifs. read( (char*)&mue_square   , sizeof(mue_square) );
        ifs. read( (char*)&sigma_square , sizeof(sigma_square) );
        ifs. read( (char*)&pos          , sizeof(pos) );

        unsigned int zpos = 0; if (pos != 0) { pos_vec.push_back(pos); pos_map[pos] = pos_vec.size(); zpos = pos_map[pos]; }

        ofs.write( (char*)&p            , sizeof(p) );
        ofs.write( (char*)&zpos         , sizeof(zpos) );
    }

    for (size_t i=0; i<n_base_faces_; ++i)
    {
        ifs.read ( (char*)&fvi[0] , sizeof(fvi[0]) );
        ifs.read ( (char*)&fvi[1] , sizeof(fvi[1]) );
        ifs.read ( (char*)&fvi[2] , sizeof(fvi[2]) );

        unsigned int fv_[3] = {fvi[0],fvi[1],fvi[2]};

        ofs.write( (char*)&fv_[0] , sizeof(fv_[0]) );
        ofs.write( (char*)&fv_[1] , sizeof(fv_[1]) );
        ofs.write( (char*)&fv_[2] , sizeof(fv_[2]) );
    }

    for (size_t i=0; i<n_details_; i++)
    {
        VsData64 vs; ifs.read( (char*)&vs , sizeof(vs) ); assert (!ifs.fail());

        if (vs.l_pos != 0) { pos_vec.push_back(vs.l_pos); pos_map[vs.l_pos] = pos_vec.size(); }
        if (vs.r_pos != 0) { pos_vec.push_back(vs.r_pos); pos_map[vs.r_pos] = pos_vec.size(); } 
    }

    for (size_t i=0; i<pos_vec.size(); i++)
    {
        const VsData64& vs = (VsData64 &) * (data + s_base_ + (pos_vec[i]-1)*s_details_);

        int zpos = 0;
        if      (vs.l_pos != 0 && vs.r_pos != 0) {zpos =  pos_map[vs.l_pos]*2+1;}
        else if (vs.l_pos != 0 && vs.r_pos == 0) {zpos =  pos_map[vs.l_pos]*2  ;}
        else if (vs.l_pos == 0 && vs.r_pos != 0) {zpos = -pos_map[vs.r_pos]*2  ;}

        unsigned long long  fund_lcut       =  vs.fund_lcut_index;
        unsigned long long  fund_rcut       =  vs.fund_rcut_index;
        unsigned int        fund_lcut_tree  =  fund_lcut >> (tree_id_bitshift_);
        unsigned int        fund_rcut_tree  =  fund_rcut >> (tree_id_bitshift_);
        unsigned int        fund_lcut_node  =  fund_lcut &  (tree_id_mask_);
        unsigned int        fund_rcut_node  =  fund_rcut &  (tree_id_mask_);

        ZipZsData64 zs; 
        zs.             p = vs.p;
        zs.fund_ccut_tree = fund_lcut_tree << 16 | fund_rcut_tree;
        zs.fund_lcut_node = fund_lcut_node;
        zs.fund_rcut_node = fund_rcut_node;
        zs.          c_pos = zpos;
        // zs.       l_pos = vs.l_pos != 0 ? pos_map[vs.l_pos] : 0;
        // zs.       r_pos = vs.r_pos != 0 ? pos_map[vs.r_pos] : 0;

        ofs.write( (char*)&zs , sizeof(zs) );
    }
    if (munmap(data, sbuf.st_size) == -1) {perror("Error un-mmapping");} close(fd);

    ifs.close(); ofs.close();

    std::cerr << n_base_vertices_ << " vertices, "
              << n_base_faces_    << " faces, "
              << n_details_       << " detail vertices\n";
}

void ZipVPM322TVPM32 (const char * ifile, const char * ofile)
{
    unsigned int        n_base_vertices_;
    unsigned int        n_base_faces_;
    unsigned int        n_details_;

    char                fileformat[16];
    unsigned int        fvi[3];
    unsigned int        pos;
    Vec3f               p, normal;
    float               radius, sin_square, mue_square, sigma_square;

    std::ifstream ifs(ifile,    std::ios::binary);
    std::ofstream ofs(ofile,    std::ios::binary); 

    if (!ifs) { std::cerr << "read error\n"; exit(1); }
    ifs.read(fileformat, 10); fileformat[10] = '\0';
    if (std::string(fileformat) != std::string("VDProgMesh"))
    { std::cerr << "Wrong file format.\n"; ifs.close(); exit(1); }

    ofs.write(fileformat, 10); // ofs << "VDProgMesh";

    ifs. read( (char*)&n_base_vertices_, sizeof(n_base_vertices_) );
    ifs. read( (char*)&n_base_faces_   , sizeof(n_base_faces_   ) );
    ifs. read( (char*)&n_details_      , sizeof(n_details_      ) );

    ofs.write( (char*)&n_base_vertices_, sizeof(n_base_vertices_) );
    ofs.write( (char*)&n_base_faces_   , sizeof(n_base_faces_   ) );
    ofs.write( (char*)&n_details_      , sizeof(n_details_      ) );

    size_t s_params_, s_base_vertices_, s_base_faces_, s_base_, s_details_;

    s_params_        = sizeof("VDProgMesh")-1+sizeof(unsigned int)+sizeof(unsigned int)+sizeof(unsigned int);
    s_base_vertices_ = sizeof(Vec3f)+sizeof(float)+sizeof(Vec3f)+sizeof(float)+sizeof(float)+sizeof(float)+sizeof(unsigned int);
    s_base_faces_    = sizeof(unsigned int)+sizeof(unsigned int)+sizeof(unsigned int);

    s_base_          = s_params_+n_base_vertices_*s_base_vertices_+n_base_faces_*s_base_faces_;
    s_details_       = sizeof(Vec3f)+sizeof(unsigned int)+sizeof(unsigned int)+sizeof(unsigned int) 
                      +sizeof(float)+sizeof(Vec3f)+sizeof(float)+sizeof(float)+sizeof(float)
                      +sizeof(float)+sizeof(Vec3f)+sizeof(float)+sizeof(float)+sizeof(float)
                      +sizeof(unsigned int)+sizeof(unsigned int);

    assert (s_details_ == sizeof(VsData32));

    // Note! We use mmap in Linux to map files and manage caches (out-of-core); (0.1.0)
    int fd; char *data; struct stat sbuf;

    if ((fd = open(ifile, O_RDONLY)) == -1)                                                     {perror("openfile");exit(1);}
    if (stat(ifile, &sbuf) == -1)                                                               {perror("filesize");exit(1);}
    if ((data = (char *) mmap(NULL, sbuf.st_size, PROT_READ, MAP_SHARED, fd, 0)) == MAP_FAILED) {perror("mmapfile");exit(1);}

    std::vector <unsigned int>  pos_vec;
    std::vector <unsigned int>  pos_map (n_base_vertices_+n_details_+1);

    for (size_t i=0; i<n_base_vertices_; ++i)
    {
        ifs. read( (char*)&p            , sizeof(p) );
        ifs. read( (char*)&radius       , sizeof(radius) );
        ifs. read( (char*)&normal       , sizeof(normal) );
        ifs. read( (char*)&sin_square   , sizeof(sin_square) );
        ifs. read( (char*)&mue_square   , sizeof(mue_square) );
        ifs. read( (char*)&sigma_square , sizeof(sigma_square) );
        ifs. read( (char*)&pos          , sizeof(pos) );

        unsigned int zpos = 0; if (pos != 0) { pos_vec.push_back(pos); pos_map[pos] = pos_vec.size(); zpos = pos_map[pos]; }

        ofs.write( (char*)&p            , sizeof(p)      );
        ofs.write( (char*)&normal       , sizeof(normal) );
        ofs.write( (char*)&zpos         , sizeof(zpos)   );
    }

    for (size_t i=0; i<n_base_faces_; ++i)
    {
        ifs.read ( (char*)&fvi[0] , sizeof(fvi[0]) );
        ifs.read ( (char*)&fvi[1] , sizeof(fvi[1]) );
        ifs.read ( (char*)&fvi[2] , sizeof(fvi[2]) );

        ofs.write( (char*)&fvi[0] , sizeof(fvi[0]) );
        ofs.write( (char*)&fvi[1] , sizeof(fvi[1]) );
        ofs.write( (char*)&fvi[2] , sizeof(fvi[2]) );
    }

    for (size_t i=0; i<n_details_; i++)
    {
        VsData32 vs; ifs.read( (char*)&vs , sizeof(vs) ); assert (!ifs.fail()); 

        if (vs.l_pos != 0) { pos_vec.push_back(vs.l_pos); pos_map[vs.l_pos] = pos_vec.size(); }
        if (vs.r_pos != 0) { pos_vec.push_back(vs.r_pos); pos_map[vs.r_pos] = pos_vec.size(); } 
    }

    for (size_t i=0; i<pos_vec.size(); i++)
    {
        const VsData32& vs = (VsData32 &) * (data + s_base_ + (pos_vec[i]-1)*s_details_);

        int zpos = 0;
        if      (vs.l_pos != 0 && vs.r_pos != 0) {zpos =  pos_map[vs.l_pos]*2+1;}
        else if (vs.l_pos != 0 && vs.r_pos == 0) {zpos =  pos_map[vs.l_pos]*2  ;}
        else if (vs.l_pos == 0 && vs.r_pos != 0) {zpos = -pos_map[vs.r_pos]*2  ;}

        ZipTVData32 zs; 
        zs.              p = vs.p;
        zs.       l_normal = vs.l_normal;
        zs.       r_normal = vs.r_normal;
        zs.fund_lcut_index = vs.fund_lcut_index;
        zs.fund_rcut_index = vs.fund_rcut_index;
        zs.          c_pos = zpos;
        // zs.       l_pos = vs.l_pos != 0 ? pos_map[vs.l_pos] : 0;
        // zs.       r_pos = vs.r_pos != 0 ? pos_map[vs.r_pos] : 0;

        ofs.write( (char*)&zs , sizeof(zs) );
    }
    if (munmap(data, sbuf.st_size) == -1) {perror("Error un-mmapping");} close(fd);

    ifs.close(); ofs.close();

    std::cerr << n_base_vertices_ << " vertices, "
              << n_base_faces_    << " faces, "
              << n_details_       << " detail vertices\n";
}

void ZipVPM642TVPM64 (const char * ifile, const char * ofile)  
{
    unsigned long long          n_base_vertices_;
    unsigned long long          n_base_faces_;
    unsigned long long          n_details_;

    char                        fileformat[16];
    unsigned long long          fvi[3];
    unsigned int                pos;
    Vec3f                       p, normal;
    float                       radius, sin_square, mue_square, sigma_square;

    std::ifstream ifs(ifile,    std::ios::binary);
    std::ofstream ofs(ofile,    std::ios::binary); 

    if (!ifs) { std::cerr << "read error\n"; exit(1); }
    ifs.read(fileformat, 10); fileformat[10] = '\0';
    if (std::string(fileformat) != std::string("VDProgMesh"))
    { std::cerr << "Wrong file format.\n"; ifs.close(); exit(1); }

    ofs.write(fileformat, 10); // ofs << "VDProgMesh";

    ifs. read( (char*)&n_base_vertices_, sizeof(n_base_vertices_) );
    ifs. read( (char*)&n_base_faces_   , sizeof(n_base_faces_   ) );
    ifs. read( (char*)&n_details_      , sizeof(n_details_      ) );

    unsigned int n_base_vertices =  n_base_vertices_;
    unsigned int n_base_faces    =  n_base_faces_;
    unsigned int n_details       =  n_details_;

    ofs.write( (char*)&n_base_vertices , sizeof(n_base_vertices ) );
    ofs.write( (char*)&n_base_faces    , sizeof(n_base_faces    ) );
    ofs.write( (char*)&n_details       , sizeof(n_details       ) );

    size_t s_params_, s_base_vertices_, s_base_faces_, s_base_, s_details_;

    s_params_        = sizeof("VDProgMesh")-1+sizeof(unsigned long long)+sizeof(unsigned long long)+sizeof(unsigned long long);
    s_base_vertices_ = sizeof(Vec3f)+sizeof(float)+sizeof(Vec3f)+sizeof(float)+sizeof(float)+sizeof(float)+sizeof(unsigned int);
    s_base_faces_    = sizeof(unsigned long long)+sizeof(unsigned long long)+sizeof(unsigned long long);

    s_base_          = s_params_+n_base_vertices_*s_base_vertices_+n_base_faces_*s_base_faces_;
    s_details_       = sizeof(Vec3f)+sizeof(unsigned int)+sizeof(unsigned long long)+sizeof(unsigned long long)+sizeof(unsigned long long) 
                      +sizeof(float)+sizeof(Vec3f)+sizeof(float)+sizeof(float)+sizeof(float)
                      +sizeof(float)+sizeof(Vec3f)+sizeof(float)+sizeof(float)+sizeof(float)
                      +sizeof(unsigned int)+sizeof(unsigned int);

    assert (s_details_ == sizeof(VsData64)); 

    // Note! We use mmap in Linux to map files and manage caches (out-of-core); (0.1.0)
    int fd; char *data; struct stat sbuf;

    if ((fd = open(ifile, O_RDONLY)) == -1)                                                     {perror("openfile");exit(1);}
    if (stat(ifile, &sbuf) == -1)                                                               {perror("filesize");exit(1);}
    if ((data = (char *) mmap(NULL, sbuf.st_size, PROT_READ, MAP_SHARED, fd, 0)) == MAP_FAILED) {perror("mmapfile");exit(1);}

    unsigned char       tree_id_bits_     = 0; while(n_base_vertices_ > ((unsigned long long) 1 << tree_id_bits_)) ++tree_id_bits_;
    unsigned char       tree_id_bitshift_ = sizeof(unsigned long long)*8   - tree_id_bits_;
    unsigned long long  tree_id_mask_     =      ((unsigned long long)-1) >> tree_id_bits_;

    std::vector <unsigned int>  pos_vec;
    std::vector <unsigned int>  pos_map (n_base_vertices_+n_details_+1);

    for (size_t i=0; i<n_base_vertices_; ++i)
    {
        ifs. read( (char*)&p            , sizeof(p) );
        ifs. read( (char*)&radius       , sizeof(radius) );
        ifs. read( (char*)&normal       , sizeof(normal) );
        ifs. read( (char*)&sin_square   , sizeof(sin_square) );
        ifs. read( (char*)&mue_square   , sizeof(mue_square) );
        ifs. read( (char*)&sigma_square , sizeof(sigma_square) );
        ifs. read( (char*)&pos          , sizeof(pos) );

        unsigned int zpos = 0; if (pos != 0) { pos_vec.push_back(pos); pos_map[pos] = pos_vec.size(); zpos = pos_map[pos]; }

        ofs.write( (char*)&p            , sizeof(p)      );
        ofs.write( (char*)&normal       , sizeof(normal) );
        ofs.write( (char*)&zpos         , sizeof(zpos)   );
    }

    for (size_t i=0; i<n_base_faces_; ++i)
    {
        ifs.read ( (char*)&fvi[0] , sizeof(fvi[0]) );
        ifs.read ( (char*)&fvi[1] , sizeof(fvi[1]) );
        ifs.read ( (char*)&fvi[2] , sizeof(fvi[2]) );

        unsigned int fv_[3] = {fvi[0],fvi[1],fvi[2]};

        ofs.write( (char*)&fv_[0] , sizeof(fv_[0]) );
        ofs.write( (char*)&fv_[1] , sizeof(fv_[1]) );
        ofs.write( (char*)&fv_[2] , sizeof(fv_[2]) );
    }

    for (size_t i=0; i<n_details_; i++)
    {
        VsData64 vs; ifs.read( (char*)&vs , sizeof(vs) ); assert (!ifs.fail());

        if (vs.l_pos != 0) { pos_vec.push_back(vs.l_pos); pos_map[vs.l_pos] = pos_vec.size(); }
        if (vs.r_pos != 0) { pos_vec.push_back(vs.r_pos); pos_map[vs.r_pos] = pos_vec.size(); } 
    }

    for (size_t i=0; i<pos_vec.size(); i++)
    {
        const VsData64& vs = (VsData64 &) * (data + s_base_ + (pos_vec[i]-1)*s_details_);

        int zpos = 0;
        if      (vs.l_pos != 0 && vs.r_pos != 0) {zpos =  pos_map[vs.l_pos]*2+1;}
        else if (vs.l_pos != 0 && vs.r_pos == 0) {zpos =  pos_map[vs.l_pos]*2  ;}
        else if (vs.l_pos == 0 && vs.r_pos != 0) {zpos = -pos_map[vs.r_pos]*2  ;}

        unsigned long long  fund_lcut       =  vs.fund_lcut_index;
        unsigned long long  fund_rcut       =  vs.fund_rcut_index;
        unsigned int        fund_lcut_tree  =  fund_lcut >> (tree_id_bitshift_);
        unsigned int        fund_rcut_tree  =  fund_rcut >> (tree_id_bitshift_);
        unsigned int        fund_lcut_node  =  fund_lcut &  (tree_id_mask_);
        unsigned int        fund_rcut_node  =  fund_rcut &  (tree_id_mask_);

        ZipTVData64 zs; 
        zs.             p = vs.p;
        zs.      l_normal = vs.l_normal;
        zs.      r_normal = vs.r_normal;
        zs.fund_ccut_tree = fund_lcut_tree << 16 | fund_rcut_tree;
        zs.fund_lcut_node = fund_lcut_node;
        zs.fund_rcut_node = fund_rcut_node;
        zs.          c_pos = zpos;
        // zs.       l_pos = vs.l_pos != 0 ? pos_map[vs.l_pos] : 0;
        // zs.       r_pos = vs.r_pos != 0 ? pos_map[vs.r_pos] : 0;

        ofs.write( (char*)&zs , sizeof(zs) );
    }
    if (munmap(data, sbuf.st_size) == -1) {perror("Error un-mmapping");} close(fd);

    ifs.close(); ofs.close();

    std::cerr << n_base_vertices_ << " vertices, "
              << n_base_faces_    << " faces, "
              << n_details_       << " detail vertices\n";
}

// Note! We read .pm files (Progressive Mesh format in OpenMesh); (0.5.3)
void PM2VPM  (const char * ifile, const char * ofile)
{
    unsigned int        n_base_vertices_;
    unsigned int        n_base_faces_;
    unsigned int        n_details_;

    char                fileformat[10];
    Vec3f               p;
    int                 v1,vl,vr,fvi[3];

    std::ifstream ifs(ifile,    std::ios::binary);
    std::ifstream ofs(ofile,    std::ios::binary);

    ifs.read(fileformat, 8); fileformat[8]='\0';
    if (std::string(fileformat) != std::string("ProgMesh")) 
    { std::cerr << "Wrong file format.\n"; ifs.close(); exit(1); }

    ifs.read( (char*)&n_base_vertices_, sizeof(n_base_vertices_) );
    ifs.read( (char*)&n_base_faces_   , sizeof(n_base_faces_   ) );
    ifs.read( (char*)&n_details_      , sizeof(n_details_      ) );

    for (size_t i=0; i<n_base_vertices_; ++i) 
    {
        ifs.read( (char*)&p , sizeof(p) ); 
    }
    for (size_t i=0; i<n_base_faces_; ++i)
    {
        ifs.read( (char*)&fvi[0] , sizeof(fvi[0]) );
        ifs.read( (char*)&fvi[1] , sizeof(fvi[1]) );
        ifs.read( (char*)&fvi[2] , sizeof(fvi[2]) );
    }
    unsigned int n_max_vertices_ = n_base_vertices_ + n_details_;
    for (size_t i=n_base_vertices_; i<n_max_vertices_; ++i) 
    {
        ifs.read( (char*)&p  , sizeof(p)  );
        ifs.read( (char*)&v1 , sizeof(v1) ); 
        ifs.read( (char*)&vl , sizeof(vl) ); 
        ifs.read( (char*)&vr , sizeof(vr) ); 
    }
    ifs.close(); ofs.close();

    std::cerr << n_base_vertices_ << " vertices, "
              << n_base_faces_    << " faces, "
              << n_details_       << " detail vertices\n";
}

bool endsWith (std::string fname, const char * suffix)
{
    return fname.rfind(suffix) != std::string::npos;
}

int main(int argc, char** argv) // Adaptive Refinement (View-dependent Rendering)
{
    // Note! We initialize all global variables here to be compatible with iOS SDK; (0.5.7) 
    if (argc == 3) { std::string ifile(argv[1]), ofile(argv[2]); 
        if      ( endsWith (ifile,"64.vpm") && endsWith (ofile, ".tvpm64") ) { ZipVPM642TVPM64 (ifile.c_str(), ofile.c_str()); } 
        else if ( endsWith (ifile,  ".vpm") && endsWith (ofile, ".tvpm"  ) ) { ZipVPM322TVPM32 (ifile.c_str(), ofile.c_str()); } 
        else if ( endsWith (ifile,"64.vpm") && endsWith (ofile, "64.zpm" ) ) { ZipVPM642ZPM64  (ifile.c_str(), ofile.c_str()); } 
        else if ( endsWith (ifile,  ".vpm") && endsWith (ofile,   ".zpm" ) ) { ZipVPM322ZPM32  (ifile.c_str(), ofile.c_str()); } 
        else if ( endsWith (ifile,   ".pm") && endsWith (ofile,   ".vpm" ) ) {         PM2VPM  (ifile.c_str(), ofile.c_str()); }
        else    { std::cerr << "Usage: ./cvmesh <FILENAME.vpm> <FILENAME.zpm> \n"; return 1; }
    }   else    { std::cerr << "Usage: ./cvmesh <FILENAME.vpm> <FILENAME.zpm> \n"; return 1; }

    return 0;

}







