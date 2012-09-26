//
//  GLTvpmMeshRenderer.cpp
//  MobiMesh
//
//  Created by Ong Yuh Shin on 5/19/11.
//  Copyright 2011 National University of Singapore. All rights reserved.
//

#include "GLTvpmMeshRenderer.h"

namespace MobiMesh {

const float GLTvpmMeshRenderer::vertex_threshold[NUM_STATES] =
    {0.01,0.001,0.0005};

GLTvpmMeshRenderer::GLTvpmMeshRenderer(GLMeshRenderer& renderer,
                                       TrulyViewDependentPM* mesh)
    : GLViewDependentMeshRenderer(renderer), mesh_(mesh),
    need_further_refine_(false), read_pixels_(true),
    read_pixels_count_(0), curr_state(NUM_STATES-1)
{
    initialize();
    set_mesh(mesh_);
    set_default_view();
    activate();
}

GLTvpmMeshRenderer::~GLTvpmMeshRenderer()
{
    deactivate();
}

void GLTvpmMeshRenderer::activate()
{
    if (mesh_)
        mesh_->activate();

    renderer_.activate();

    // Reload the color buffer after activating the renderer_
    std::vector<OpenMesh::Vec4uc> mcolors_(
        TrulyViewDependentPM::MAX_BUFF_FACES*3);
    unsigned char r = 0, g = 0, b = 1;

    // Assign each face (up to MAX_BUFF_FACES) a unique color.
    // We will use this color to determine how many faces are on
    // screen and how much of the face can be seen.
    // The mapping of face color to face index is made 1-1 so
    // that we can retrieve the corresponding face index from the
    // color later in the adapt_mesh() method.
    for (size_t i = 0; i < TrulyViewDependentPM::MAX_BUFF_FACES; i++)  
    {
        mcolors_[i*3  ] = OpenMesh::Vec4uc(r,g,b,255);
        mcolors_[i*3+1] = OpenMesh::Vec4uc(r,g,b,255); 
        mcolors_[i*3+2] = OpenMesh::Vec4uc(r,g,b,255);
		
        if (b == 255)
        {
            b = 0;
            if (g == 255)
            {
                g = 0;
#ifdef DEBUG
                if (r == 255)
                    std::cerr<<"overflow, insufficient colors to make all "
                    "faces unique"<<std::endl;
#endif
                r++;
            }
            else
                g++;
        }
        else
            b++;
    }

    renderer_.set_color(mcolors_);
    image_ = new unsigned char[get_viewport_width() *
                               get_viewport_height() * 4];
}

void GLTvpmMeshRenderer::deactivate()
{
    renderer_.deactivate();

    if (mesh_)
        mesh_->deactivate();

    if (image_)
    {
        delete image_;
        image_ = NULL;
    }
}

void GLTvpmMeshRenderer::initialize()
{
    ec_threshold_ = vertex_threshold[NUM_STATES-1] *
        get_viewport_width() * get_viewport_height();
    ec_threshold_ = (sqrt(ec_threshold_) / 2) - 1;
    ec_threshold_ = (ec_threshold_ < 1) ? 1 : ec_threshold_;
}

void GLTvpmMeshRenderer::read_pixels()
{
    // Render each face with a specific color first
    renderer_.render_with_color(true);
    renderer_.render();

    // Read the pixels from the rendered screen back so that we can
    // determine which colors (hence faces) were rendered, how big
    // their surface area is, and then decide on the vertex weights for
    // the mesh refinement
    glReadPixels(0,0, get_viewport_width(), get_viewport_height(),
                GL_RGBA, GL_UNSIGNED_BYTE, image_);
    //glReadPixels(0, 0, get_viewport_width(), get_viewport_height(),
    //             GL_IMPLEMENTATION_COLOR_READ_FORMAT,GL_IMPLEMENTATION_COLOR_READ_TYPE, image_);
    
    renderer_.render_with_color(false);

    ++read_pixels_count_;
}

void GLTvpmMeshRenderer::update_vertex_weights()
{
    mesh_->reset_vertex_weights();
    
    // Read the pixels from the image array and update the vertex weights
    // accordingly. First find the weight of each face based on the number of
    // pixels containing that face color.
    size_t num_faces = mesh_->num_faces();
    std::vector<unsigned int> mweight(num_faces, 0);
    for (int i = 0; i < get_viewport_width() * get_viewport_height(); i++) 
    {
        unsigned char  r = image_[i*4  ];
        unsigned char  g = image_[i*4+1];
        unsigned char  b = image_[i*4+2];
        unsigned int seq = r * 65536 + g * 256 + b;

#ifdef DEBUG
        assert(seq <= num_faces);
#endif

        if (seq != 0)
            mweight[seq-1]++;
    }

    // Assign the weights of the faces to the corresponding vertices making up
    // the face
    for (size_t i = 0; i < mweight.size(); i++)
    {
        unsigned int weight = mweight[i]/3;
        if (weight != 0)
        {
            mesh_->add_vertex_weight(mesh_->face_index_to_vertex_index(3*i),
                                     weight);
            mesh_->add_vertex_weight(mesh_->face_index_to_vertex_index(3*i+1),
                                     weight);
            mesh_->add_vertex_weight(mesh_->face_index_to_vertex_index(3*i+2),
                                     weight);
        }
    }
}

bool GLTvpmMeshRenderer::adapt_mesh()
{
    if (!mesh_)
        return false; // no mesh to render

    if (view_changed())
    {
        // Everytime we are at a new viewpoint, we reset the flags and begin by
        // first reading the pixels to determine the weight of each face.
        read_pixels_ = true;
        need_further_refine_ = false;
        read_pixels_count_ = 0;

        // We cycle around the possible states, lowering the vs threshold each
        // time until we reach the minimum, then we go repeat the cycle and
        // start at the maximum.
        ++curr_state;
        curr_state = curr_state % NUM_STATES;
        need_further_refine_ = true;
    }
    else if (curr_state != NUM_STATES - 1)
    {
        // Continue to the next state until we reach VS_THRESHOLD_4. Here,
        // instead of looping back to VS_THRESHOLD_100 as we do when
        // view_changed() is true, we stay at a low threshold since the
        // current view is likely to have very fine faces so we need a low
        // threshold to refine even more detail.
        ++curr_state;
        need_further_refine_ = true;
    }

    bool is_adapted = false;

    // We must call read_pixels() before any refinement is made.
    // If we perform the refinement first, then read pixels within the same
    // frame, the mesh will be updated in the mesh object. However, the VBOs in
    // renderer_ will not be updated and reading back the pixels, which require
    // rendering, and will render the old set of vertices from before
    // refinement.
    if (read_pixels_)
    {
        read_pixels();
        update_vertex_weights();
        read_pixels_ = false;
        need_further_refine_ = true;
    }

    if (view_changed() || need_further_refine_)
    {
        int vs_threshold = vertex_threshold[curr_state] *
            get_viewport_width() * get_viewport_height();
        vs_threshold = sqrt(vs_threshold);
        vs_threshold = (vs_threshold < 3) ? 3 : vs_threshold;

        need_further_refine_ = mesh_->adapt(MAX_REFINEMENTS,
                                            ec_threshold_, vs_threshold);
        is_adapted = true;
    }

    // Once mesh refinement has converged for the pixels we read back, we will
    // read_pixels_ again to determine the faces that are on screen based
    // on the converged mesh. This is because one read of the pixels do
    // not give an accurate enough estimation of the face/vertex weights
    // as the initial mesh is very coarse.
    // However, we limit also the number of times we read the pixels back
    // for each viewpoint so that the mesh refinement can converge.
    if (!need_further_refine_ && read_pixels_count_ < 3)
        read_pixels_ = true;

    return is_adapted;
}

}