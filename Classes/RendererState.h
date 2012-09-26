//
//  RendererState.h
//  MobiMesh
//
//  Created by Ong Yuh Shin on 5/26/11.
//  Copyright 2011 National University of Singapore. All rights reserved.
//

#import <Foundation/Foundation.h>

#import "glm.hpp"
#import "GLMeshRendererFactory.h"


@interface RendererState : NSObject {

    // Factory for creating the renderer
    MobiMesh::GLMeshRendererFactory* rendererFactory;
    // Handle to the renderer object created by factory
    MobiMesh::GLMeshRendererFactory::RendererHandle rendererHandle;
    // Mesh filename for loading/reloading the mesh
    NSString* meshFilename;
    NSString* vertexShader;
    NSString* fragmentShader;

    // The last lookat/transformation matrices before the
    // renderer was unloaded. To be applied once renderer is reloaded.
    glm::mat4 lookatMatrix;
    glm::mat4 transformationMatrix;

    // Loading mesh for first time?
    bool firstTimeLoading;
}

@property (readonly, nonatomic) MobiMesh::GLMeshRendererFactory* rendererFactory;
@property (readonly, nonatomic) MobiMesh::GLMeshRendererFactory::RendererHandle rendererHandle;
@property (readonly, nonatomic) NSString* meshFilename;
@property (readonly, nonatomic) NSString* vertexShader;
@property (readonly, nonatomic) NSString* fragmentShader;

// Creates the renderer object using the given factory.
// Mesh/renderer is not loaded upon initialization. Use loadRenderer for that.
// Use this if using OpenGLES1
- (id)initWithFactory: (MobiMesh::GLMeshRendererFactory*) factory
         withFilename: (NSString*) filename;

// Same as initWithFactory: ... withFilename: ... but for OpenGL ES2
// The vertex and fragment shader must be supplied and cannot be nil
- (id)initWithFactory: (MobiMesh::GLMeshRendererFactory*) factory
         withFilename: (NSString*) filename
     withVertexShader: (NSString*) vshader
   withFragmentShader: (NSString*) fshader;

// Loads the renderer. No effect if renderer is already loaded.
// Throws UnableToCreateRendererException or InvalidMeshTypeException
// (propagated from GLMeshRendererFactory)
- (void)loadRendererWithViewportWidth: (GLsizei) width
                   withViewportHeight:(GLsizei) height;

// Unloads the renderer. No effect if renderer is already loaded.
- (void)unloadRenderer;

// Whether the renderer is currently loaded in memory
- (BOOL)isLoaded;

@end
