//
//  RendererState.m
//  MobiMesh
//
//  Created by Ong Yuh Shin on 5/26/11.
//  Copyright 2011 National University of Singapore. All rights reserved.
//

#import "RendererState.h"


@implementation RendererState

@synthesize rendererFactory, rendererHandle, meshFilename, vertexShader, fragmentShader;

- (id)initWithFactory:(MobiMesh::GLMeshRendererFactory *)factory
         withFilename:(NSString *)filename
{
    self = [super init];
    if (self)
    {
        rendererFactory = factory;
        rendererHandle = MobiMesh::GLMeshRendererFactory::INVALID_HANDLE;
        meshFilename = [filename retain];
        vertexShader = nil;
        fragmentShader = nil;
        firstTimeLoading = true;
    }
    return self;
}

- (id)initWithFactory:(MobiMesh::GLMeshRendererFactory *)factory
         withFilename:(NSString *)filename
     withVertexShader:(NSString *)vshader
   withFragmentShader:(NSString *)fshader
{
    self = [super init];
    if (self)
    {
        rendererFactory = factory;
        rendererHandle = MobiMesh::GLMeshRendererFactory::INVALID_HANDLE;
        meshFilename = [filename retain];
        vertexShader = [vshader retain];
        fragmentShader = [fshader retain];
        firstTimeLoading = true;
    }
    return self;
}

- (void)dealloc
{
    [meshFilename release];
    if (vertexShader != nil)
        [vertexShader release];
    if (fragmentShader != nil)
        [fragmentShader release];
    [super dealloc];
}

- (void)loadRendererWithViewportWidth:(GLsizei)width
                   withViewportHeight:(GLsizei)height
{
    // Renderer is active
    if ([self isLoaded])
        return;

    if (vertexShader == nil) // using ES1
        rendererHandle = rendererFactory->new_mesh_renderer(
            [meshFilename cStringUsingEncoding:NSASCIIStringEncoding],
            width, height);
    else
    {
        // Using ES2
        rendererHandle =
            rendererFactory->new_es2_mesh_renderer(
                [meshFilename cStringUsingEncoding:NSASCIIStringEncoding],
                [vertexShader cStringUsingEncoding:NSUTF8StringEncoding],
                [fragmentShader cStringUsingEncoding:NSUTF8StringEncoding],
                width, height);
    }

    if (!firstTimeLoading)
    {
        MobiMesh::IMeshRenderer* renderer = rendererFactory->get_renderer(rendererHandle);
        if (renderer) // load the saved view
            renderer->set_view(lookatMatrix, transformationMatrix);
    }

    firstTimeLoading = false;
}

- (void)unloadRenderer
{
    // Save the view in case we need to reload the renderer/mesh later
    MobiMesh::IMeshRenderer* renderer = rendererFactory->get_renderer(rendererHandle);
    if (renderer)
        renderer->capture_view(lookatMatrix, transformationMatrix);

    // Free up the renderer and set the handle back to invalid.
    rendererFactory->free_mesh_renderer(rendererHandle);
    rendererHandle = MobiMesh::GLMeshRendererFactory::INVALID_HANDLE;
}

- (BOOL)isLoaded
{
    return (rendererFactory->get_renderer(rendererHandle) != NULL);
}

@end
