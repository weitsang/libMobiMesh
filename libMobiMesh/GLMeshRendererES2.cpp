//
//  GLMeshRendererES2.cpp
//  MobiMesh
//
//  Created by Ong Yuh Shin on 4/9/11.
//  Copyright 2011 National University of Singapore. All rights reserved.
//

#include "GLMeshRendererES2.h"
#include "InvalidGLProgramException.h"
#include "ProgramLinkingException.h"
#include "ShaderCompileException.h"

// glm::vec3, glm::vec4, glm::ivec4, glm::mat4
#include "glm.hpp"
// glm::translate, glm::rotate, glm::scale
#include "matrix_transform.hpp"

namespace MobiMesh {

const char* GLMeshRendererES2::UNIFORM_MODELVIEW_PROJECTION_STRING = "modelviewProjection";
const char* GLMeshRendererES2::UNIFORM_NORMAL_MATRIX_STRING = "normalMatrix";
const char* GLMeshRendererES2::UNIFORM_RENDER_WITH_COLOR_STRING = "renderWithColor";

const char* GLMeshRendererES2::ATTRIBUTE_POSITION_STRING = "position";
const char* GLMeshRendererES2::ATTRIBUTE_NORMAL_STRING = "normal";
const char* GLMeshRendererES2::ATTRIBUTE_COLOR_STRING = "color";

GLMeshRendererES2::GLMeshRendererES2(
    const char* vertex_shader,
    const char* fragment_shader,
    GLsizei     viewport_width,
    GLsizei     viewport_height)
    : program_(0),
    GLMeshRenderer(viewport_width, viewport_height)
{
    initialize(vertex_shader, fragment_shader);
}

GLMeshRendererES2::GLMeshRendererES2(
    const char* vertex_shader,
    const char* fragment_shader,
    IMesh*      mesh,
    GLsizei     viewport_width,
    GLsizei     viewport_height)
	: program_(0),
    GLMeshRenderer(mesh, viewport_width, viewport_height)
{
    initialize(vertex_shader, fragment_shader);
}

GLMeshRendererES2::~GLMeshRendererES2()
{
    if (program_)
        glDeleteProgram(program_);
}

void GLMeshRendererES2::initialize(
    const char *vertex_shader, const char *fragment_shader)
{
    program_ = build_program(vertex_shader, fragment_shader);

    // Commented out because it causes EXC_BAD_EXCEPTION to be thrown
    // for some unknown reason and random times when a new renderer
    // is created.
    // validate_program(program_);

    // Initialize uniform variables accordingly
    // (matrices will take a default identity matrix)
    glm::mat4 id(1.0f);

    glUseProgram(program_);

    uniform_modelview_projection_ = glGetUniformLocation(program_, UNIFORM_MODELVIEW_PROJECTION_STRING);
    glUniformMatrix4fv(uniform_modelview_projection_, 1, 0, &id[0][0]);

    uniform_normal_matrix_ = glGetUniformLocation(program_, UNIFORM_NORMAL_MATRIX_STRING);
    glUniformMatrix4fv(uniform_normal_matrix_, 1, 0, &id[0][0]);

    uniform_color_ = glGetUniformLocation(program_, UNIFORM_RENDER_WITH_COLOR_STRING);
    glUniform1i(uniform_color_, 0);
}

void GLMeshRendererES2::setup_modelview_projection()
{
    // compute projection/perspective matrix
    glm::mat4 perspective_matrix =
        glm::perspective((GLfloat) GLMeshRenderer::FOVY,
                         (GLfloat) get_viewport_width() / (GLfloat) get_viewport_height(),
                         view_params_.get_near_plane(),
                         view_params_.get_far_plane());

    // compute modelview matrix
    glm::mat4 lookat, transformation;
    capture_view(lookat, transformation);
    glm::mat4 modelview_matrix = lookat * transformation;

    // multiply the 2 matrices for the modelview_projection (note: projection before modelview)
    glm::mat4 modelview_projection_matrix = perspective_matrix * modelview_matrix;

    glUniformMatrix4fv(uniform_modelview_projection_, 1, 0, &modelview_projection_matrix[0][0]);

    // Also setup the normal matrix (since this is affected by the modelview matrix)
    glm::mat4 normal_matrix = glm::transpose(glm::inverse(modelview_matrix));
    glUniformMatrix4fv(uniform_normal_matrix_, 1, 0, &normal_matrix[0][0]);
}

void GLMeshRendererES2::setup_vertices(const GLvoid* pointer)
{
    glVertexAttribPointer(ATTRIBUTE_POSITION, 3, GL_FLOAT, GL_FALSE, 0, pointer);
}

void GLMeshRendererES2::setup_normals(const GLvoid* pointer)
{
    glVertexAttribPointer(ATTRIBUTE_NORMAL, 3, GL_FLOAT, GL_TRUE, 0, pointer);
}

void GLMeshRendererES2::setup_colors(const GLvoid *pointer)
{
    glVertexAttribPointer(ATTRIBUTE_COLOR, 4, GL_UNSIGNED_BYTE, GL_TRUE, 0, pointer);

#ifdef DEBUG
    validate_program(program_);
#endif
}

void GLMeshRendererES2::rendering_setup()
{
    glEnableVertexAttribArray(ATTRIBUTE_POSITION);

    // set the uniform color flag accordingly
    if (render_with_color())
    {
        glEnableVertexAttribArray(ATTRIBUTE_COLOR);
        glUniform1i(uniform_color_, 1);
    }
    else
    {
        glEnableVertexAttribArray(ATTRIBUTE_NORMAL);
        glUniform1i(uniform_color_, 0);
    }
}

void GLMeshRendererES2::rendering_cleanup()
{
    glDisableVertexAttribArray(ATTRIBUTE_POSITION);

    if (render_with_color())
        glDisableVertexAttribArray(ATTRIBUTE_COLOR);
    else
        glDisableVertexAttribArray(ATTRIBUTE_NORMAL);
}

GLuint GLMeshRendererES2::build_program(
    const char *vertex_shader_src,
    const char *fragment_shader_src) const
{
    // build shader might throw ShaderCompileException, let exception propagate
    GLuint vertex_shader = build_shader(vertex_shader_src, GL_VERTEX_SHADER);
    GLuint fragment_shader = build_shader(fragment_shader_src, GL_FRAGMENT_SHADER);

    GLuint program_handle = glCreateProgram();
    if (!program_handle)
    {
        // Unable to create program, free shaders and throw exception
        if (vertex_shader)
            glDeleteShader(vertex_shader);
        if (fragment_shader)
            glDeleteShader(fragment_shader);
        throw InvalidGLProgramException("Unable to create program");
    }

    // Attach shaders
    glAttachShader(program_handle, vertex_shader);
    glAttachShader(program_handle, fragment_shader);

    // Bind attribute locations
    glBindAttribLocation(program_handle, ATTRIBUTE_POSITION, ATTRIBUTE_POSITION_STRING);
    glBindAttribLocation(program_handle, ATTRIBUTE_NORMAL, ATTRIBUTE_NORMAL_STRING);
    glBindAttribLocation(program_handle, ATTRIBUTE_COLOR, ATTRIBUTE_COLOR_STRING);

    glLinkProgram(program_handle);

    // Free shaders
    if (vertex_shader)
        glDeleteShader(vertex_shader);
    if (fragment_shader)
        glDeleteShader(fragment_shader);

    GLint success;
    glGetProgramiv(program_handle, GL_LINK_STATUS, &success);
    if (success == GL_FALSE)
    {
        GLchar messages[256];
        glGetProgramInfoLog(program_handle, sizeof(messages), 0, &messages[0]);
        if (program_handle)
            glDeleteProgram(program_handle);

        throw ProgramLinkingException(messages);
    }

    return program_handle;
}

void GLMeshRendererES2::validate_program(GLuint program) const
{
    glValidateProgram(program);

    GLint success;
    glGetProgramiv(program, GL_VALIDATE_STATUS, &success);
    if (success == GL_FALSE)
    {
        GLchar messages[256];
        glGetProgramInfoLog(program, sizeof(messages), 0, &messages[0]);
        if (program_)
            glDeleteProgram(program_);
#ifdef DEBUG
        std::cout << "Unable to validate program, reason: "
                  << messages << std::endl;
#endif
        throw InvalidGLProgramException(messages);
    }
}

GLuint GLMeshRendererES2::build_shader(const char *src, GLenum shader_type) const
{
    GLuint shader_handle = glCreateShader(shader_type);
    glShaderSource(shader_handle, 1, &src, 0);
    glCompileShader(shader_handle);

    GLint success;
    glGetShaderiv(shader_handle, GL_COMPILE_STATUS, &success);
    
    if (success == GL_FALSE)
    {
        GLchar messages[256];
        glGetShaderInfoLog(shader_handle, sizeof(messages), 0, &messages[0]);
        glDeleteShader(shader_handle);
        throw ShaderCompileException(src, messages);
    }

    return shader_handle;
}

}