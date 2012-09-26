//
//  MobiMeshViewController.h
//  MobiMesh
//
//  Created by Ong Yuh Shin on 3/4/11.
//  Copyright 2011 National University of Singapore. All rights reserved.
//

#import <UIKit/UIKit.h>

#import <OpenGLES/EAGL.h>

#import <OpenGLES/ES1/gl.h>
#import <OpenGLES/ES1/glext.h>
#import <OpenGLES/ES2/gl.h>
#import <OpenGLES/ES2/glext.h>

#import <list>

#import "glm.hpp"

#import "MenuViewController.h"
#import "ModelViewCaptureContainer.h"
#import "RendererState.h"

#import "GLMeshRendererFactory.h"
#import "GLMeshRenderer.h"
#import "TrackBall.h"

@interface MobiMeshViewController : UIViewController
{
    EAGLContext *context;
    GLuint program;

    BOOL animating;
    NSInteger animationFrameInterval;
    CADisplayLink *displayLink;

	// UI objects
    IBOutlet UILabel* meshFilename;
    IBOutlet UILabel* errorMessageLabel;
    IBOutlet UIActivityIndicatorView* meshLoadingIndicator;

    MobiMesh::GLMeshRendererFactory rendererFactory;

    // For switching between meshes/renderers
    // Meshes will ordered in (pm, vpm, tvpm) triplets
    // i.e. The PM format of mesh i is found at meshNames[3i]
    // and the VPM format is at meshNames[3i+1], and TVPM at meshNames[3i+2]
    NSMutableArray* meshNames; // mesh names to be loaded

    // Renderers[i] refers to the mesh renderer of mesh at meshNames[i]
    // Each is a RendererState object.
    std::vector<RendererState*> renderers;
    int currRendererIndex;

    // For capturing specific views for each mesh.
    // Capturing a view for mesh i will capture that view for the pm,vpm,tvpm
    // versions.
    // The view is stored as a pair of matrices, first is the lookat matrix,
    // second is the transformation matrix.
    // When a mesh is loaded, the default view upon loading the mesh object
    // is captured if there are no existing captures.
    // capturedView[i] corresponds to a view for meshNames[3i + [0,2]]
    ModelViewCaptureContainer* capturedView;
    int numFormats;

    // Flag to determine whether a mesh is loading. All UI methods must not
    // proceed if this mesh is loading (i.e. if flag is set, the methods
    // simply return without doing anything) This is to prevent there being
    // multiple threads loading renderers at the same time. The renderer
    // loading is not thread safe.
    bool isMeshLoading;

    // Whether to autorotate the mesh (rotation takes place about y axis)
    bool autoRotate;

	// Touch-handling
	float      startingTouchDistance, previousScale;
	CGPoint    lastMovementPosition, previousDirectionOfPanning;
	BOOL       twoFingersAreMoving, pinchGestureUnderway;

    // Shaders when using ES2
    NSString* vertexShaderSource;
    NSString* fragmentShaderSource;

	///// MobiMesh: end /////
}

@property (readonly, nonatomic, getter=isAnimating) BOOL animating;
@property (nonatomic) NSInteger animationFrameInterval;

@property (retain, readonly, nonatomic) NSArray* meshNames;
@property (retain, nonatomic) ModelViewCaptureContainer* capturedView;

@property (retain, nonatomic) UILabel* meshFilename;
@property (retain, nonatomic) UILabel* errorMessageLabel;
@property (retain, nonatomic) UIActivityIndicatorView* meshLoadingIndicator;

- (void)startAnimation;
- (void)stopAnimation;

// Selects the active mesh to be rendered and loads the mesh if necessary
// Modifies currRendererIndex to be the supplied rendererIndex
- (void)selectMeshAtIndex:(int)rendererIndex;
// Unloads the mesh/renderer at the given index
- (void)unloadMeshAtIndex:(int)rendererIndex;

// Gets the current mesh renderer, NULL if the renderer has not been
// created.
- (MobiMesh::IMeshRenderer*)getCurrentRenderer;
- (MobiMesh::IMeshRenderer*)getRendererAtIndex:(int)rendererIndex;

// For saving/loading captured viewpoints
- (NSString*)dataFilePath;
// This will save the captured viewpoints when the application quits
- (void)saveCapturedData;

// UI actions
- (IBAction)showMenu:(id)sender;
- (IBAction)toggleAutoRotate:(id)sender;
- (IBAction)captureView:(id)sender;
- (IBAction)loadView:(id)sender;
- (IBAction)previousMesh:(id)sender;
- (IBAction)nextMesh:(id)sender;

// Loads the mesh and creates the renderer accordingly
// This runs on the main thread and is made to do so by loadMeshInBackground
- (void)loadMeshAtIndex:(int)index;
// Sets error string. Since it modifies a UILabel object, it should be
// called on the main thread.
- (void)setErrorString:(NSString*)error;

// Actions to take before the mesh starts loading (must run on main thread)
- (void)willStartMeshLoading;
// Actions to take once the mesh stops loading (must run on main thread)
- (void)didStopMeshLoading;

// Captures the view for the current mesh and adds it to the list of
// possible views for this mesh
- (void)captureCurrentView;
- (void)captureViewAtRendererIndex:(int)rendererIndex;
// Loads the view at viewIndex for the current mesh
- (void)loadViewAtViewIndex:(int)viewIndex;

// A background thread to load the mesh and set the loading indicator
// accordingly
- (void)loadMeshInBackground;

// Calls loadMeshInBackground and lets it run in the background if the current
// mesh is not loaded
- (void)loadCurrMeshIfNotLoaded;

@end
