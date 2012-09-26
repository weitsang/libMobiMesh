//
//  TVPMRendererViewController.h
//  MobiMesh
//
//  Created by Ong Yuh Shin on 6/30/11.
//  Copyright 2011 National University of Singapore. All rights reserved.
//

#import <UIKit/UIKit.h>

#import <OpenGLES/EAGL.h>

#import <OpenGLES/ES1/gl.h>
#import <OpenGLES/ES1/glext.h>
#import <OpenGLES/ES2/gl.h>
#import <OpenGLES/ES2/glext.h>

#import "GLMeshRenderer.h"
#import "GLTvpmMeshRenderer.h"
#import "TrulyViewDependentPM.h"

@interface TVPMRendererViewController : UIViewController {
@private
    EAGLContext *context;
    GLuint program;
    
    BOOL animating;
    NSInteger animationFrameInterval;
    CADisplayLink *displayLink;

    // *** MobiMesh Tutorial: Declare necessary renderer objects
    MobiMesh::GLMeshRenderer* renderer;
    MobiMesh::GLTvpmMeshRenderer* tvpmRenderer;
    MobiMesh::TrulyViewDependentPM* mesh;
}

@property (readonly, nonatomic, getter=isAnimating) BOOL animating;
@property (nonatomic) NSInteger animationFrameInterval;

- (void)startAnimation;
- (void)stopAnimation;

@end
