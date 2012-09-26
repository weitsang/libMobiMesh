//
//  GLMeshRendererES2.h
//  MobiMesh
//
//  Created by Ong Yuh Shin on 4/9/11.
//  Copyright 2011 National University of Singapore. All rights reserved.
//

#ifndef __GL_MESH_RENDERER_ES2_H_
#define __GL_MESH_RENDERER_ES2_H_

#include <OpenGLES/ES2/gl.h>
#include <OpenGLES/ES2/glext.h>

#include "glm.hpp"

#include "IMesh.h"
#include "IMeshRenderer.h"
#include "GLMeshRenderer.h"
#include "GLMeshRendererViewingParameters.h"

namespace MobiMesh {

/**
 * Concrete class responsible for making OpenGL ES 2.0 specific calls
 * when rendering.
 */
class GLMeshRendererES2 : public GLMeshRenderer 
{
public:

	///////////// Constructors

    /**
     * Creates a renderer with the supplied shaders. The renderer uses
     * OpenGL ES 2.0 for rendering.
     *
     * @throw ShaderCompileException if shaders could not be built.
     * @throw ProgramLinkingException if program could not be linked.
     * @throw InvalidGLProgramException if program could not be validated.
     * @param[in] vertex_shader string containing the vertex shader code
     * (not the shader's filename)
     * @param[in] fragment_shader string containing the fragment shader code
     * (not the shader's filename)
     * @param[in] viewport_width the viewport's width
     * @param[in] viewport_height the viewport's height
     */
	GLMeshRendererES2(
        const char* vertex_shader,
        const char* fragment_shader,
        GLsizei viewport_width,
        GLsizei viewport_height);

    /**
     * Creates a renderer with the supplied shaders and the supplied
     * mesh to be rendered. The renderer uses OpenGL ES 2.0 for rendering.
     *
     * @throw ShaderCompileException if shaders could not be built.
     * @throw ProgramLinkingException if program could not be linked.
     * @throw InvalidGLProgramException if program could not be validated.
     * @param[in] vertex_shader string containing the vertex shader code
     * (not the shader's filename)
     * @param[in] fragment_shader string containing the fragment shader code
     * (not the shader's filename)
     * @param[in] mesh mesh object to render. Caller is responsible for
     * ensuring that lifetime of mesh object exceeds that of the created
     * renderer.
     * @param[in] viewport_width the viewport's width
     * @param[in] viewport_height the viewport's height
     */
	GLMeshRendererES2(
        const char* vertex_shader,
        const char* fragment_shader,
        IMesh*      mesh,
        GLsizei     viewport_width,
        GLsizei     viewport_height);

	virtual ~GLMeshRendererES2();

    ///////////// Methods inherited from GLMeshRenderer

    void render() { // glUseProgram(program_); 
        GLMeshRenderer::render(); }

    ///////////// Methods specific to OpenGL ES 2.0 rendering

    GLuint get_program() { return program_; }

    ///////////// Uniform/attribute variable names for shaders

    static const char* UNIFORM_MODELVIEW_PROJECTION_STRING;
    static const char* UNIFORM_NORMAL_MATRIX_STRING;
    static const char* UNIFORM_RENDER_WITH_COLOR_STRING;

    // Attributes for each vertex's position and normal in object space
    static const char* ATTRIBUTE_POSITION_STRING;
    static const char* ATTRIBUTE_NORMAL_STRING;
    static const char* ATTRIBUTE_COLOR_STRING;

    ///////////// Attribute index location for shaders

    enum attribute_index
    {
        ATTRIBUTE_POSITION,
        ATTRIBUTE_NORMAL,
        ATTRIBUTE_COLOR,
        NUM_ATTRIBUTES
    };

protected:

    ///////////// Methods inherited from GLMeshRenderer
    
    virtual void setup_modelview_projection();
    virtual void setup_vertices(const GLvoid* pointer);
    virtual void setup_normals(const GLvoid* pointer);
    virtual void setup_colors(const GLvoid* pointer);
    virtual void rendering_setup();
    virtual void rendering_cleanup();

private:
    
	GLMeshRendererES2(const GLMeshRendererES2&);
	GLMeshRendererES2& operator=(const GLMeshRendererES2&);
    
    /**
     * Called by the constructor to initialize the program (build and link
     * shaders).
     * @throw ShaderCompileException if shaders could not be built.
     * @throw ProgramLinkingException if program could not be linked.
     * @throw InvalidGLProgramException if program could not be validated.
     */
    void initialize(const char* vertex_shader, const char* fragment_shader);

    /**
     * Builds a program comprising of the vertex and fragment shaders.
     * @param[in] vertex_shader_src src for the vertex shader.
     * @param[in] fragment_shader_src src for the fragment shader.
     * @return program handle of the linked program
     * @throw ShaderCompileException if shaders could not be built.
     * @throw ProgramLinkingException if program could not be linked.
     */
    GLuint build_program(
        const char* vertex_shader_src,
        const char* fragment_shader_src) const;
    
    /**
     * Validates the given program.
     * @param[in] program handle of program to be validated
     * @throw InvalidGLProgramException if program is invalid.
     */
    void validate_program(GLuint program) const;

    /**
     * Builds the given shader.
     * @param[in] src src for the shader file
     * @param[in] the shader type to be supplied to the glCreateShader call
     * @return shader handle of the compiled shader
     * @throw ShaderCompileException if unable to build shader
     */
    GLuint build_shader(const char* src, GLenum shader_type) const;

    /// Linked shaders for OpenGL ES 2.0
    GLuint program_;

    /// Uniform locations
    GLint uniform_modelview_projection_, uniform_normal_matrix_, uniform_color_;
};

}

#endif
