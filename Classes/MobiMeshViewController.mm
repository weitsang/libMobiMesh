//
//  MobiMeshViewController.m
//  MobiMesh
//
//  Created by Ong Yuh Shin on 3/4/11.
//  Copyright 2011 National University of Singapore. All rights reserved.
//

#import <QuartzCore/QuartzCore.h>

#import "MobiMeshAppDelegate.h"
#import "MobiMeshViewController.h"
#import "EAGLView.h"

#import "GLMeshRenderer.h"
#import "GLMeshRendererFactory.h"
#import "ViewDependentProgressiveMesh.h"
#import "InvalidMeshTypeException.h"
#import "UnableToCreateRendererException.h"

#define kViewCaptureFilename @"views"
#define kSavedView @"savedView"

@interface MobiMeshViewController ()
@property (nonatomic, retain) EAGLContext *context;
@property (nonatomic, assign) CADisplayLink *displayLink;
@end

@implementation MobiMeshViewController

@synthesize animating, context, displayLink;
@synthesize meshNames, capturedView;
@synthesize meshFilename;
@synthesize errorMessageLabel;
@synthesize meshLoadingIndicator;

- (void)awakeFromNib
{
    EAGLContext *aContext = [[EAGLContext alloc] initWithAPI:kEAGLRenderingAPIOpenGLES2];

    if (!aContext) // Fall back to ES1
        aContext = [[EAGLContext alloc] initWithAPI:kEAGLRenderingAPIOpenGLES1];

    if (!aContext)
        NSLog(@"Failed to create ES context");
    else if (![EAGLContext setCurrentContext:aContext])
        NSLog(@"Failed to set ES context current");
    
	self.context = aContext;
	[aContext release];

	[(EAGLView *)self.view setContext:context];
	[(EAGLView *)self.view setFramebuffer];

	if ([context API] == kEAGLRenderingAPIOpenGLES2)
    {
        // Load the shaders if using ES2
        NSString *vertexShaderFile = [[NSBundle mainBundle] pathForResource:@"Shader" ofType:@"vsh"];
        NSString *fragmentShaderFile = [[NSBundle mainBundle] pathForResource:@"Shader" ofType:@"fsh"];
        
        vertexShaderSource =
            [[NSString stringWithContentsOfFile:vertexShaderFile
              encoding:NSUTF8StringEncoding error:nil] retain];
        if (!vertexShaderSource)
        {
            NSLog(@"Failed to load vertex shader");
            exit(0);
        }

        fragmentShaderSource =
            [[NSString stringWithContentsOfFile:fragmentShaderFile
            encoding:NSUTF8StringEncoding error:nil] retain];
        if (!fragmentShaderSource)
        {
            NSLog(@"Failed to load fragment shader");
            exit(0);
        }
    }
    else
    {
        vertexShaderSource = nil;
        fragmentShaderSource = nil;
    }
    
	animating = FALSE;
	animationFrameInterval = 1;
	self.displayLink = nil;

    // Init the names of meshes to be viewable
    meshNames = [[NSArray alloc] initWithObjects:
                 @"thai.pm", @"thai.vpm", @"thai.tvpm",
                 @"horse.pm", @"horse.vpm", @"horse.tvpm64",
                 @"happy.pm", @"happy.vpm", @"happy.tvpm",
                 @"lucy.pm", @"lucy.vpm", @"lucy.tvpm64",
                 nil];
    numFormats = 3;

    // Init vector sizes
    int numMeshes = [meshNames count];
    renderers.resize(numMeshes, nil);

    // Set mesh to auto rotate initially
    autoRotate = true;

    // Load the captured modelview matrices if any
    NSString* filePath = [self dataFilePath];
    if ([[NSFileManager defaultManager] fileExistsAtPath:filePath])
    {
        NSData* data = [[NSMutableData alloc] initWithContentsOfFile:[self dataFilePath]];
        NSKeyedUnarchiver* unarchiver = [[NSKeyedUnarchiver alloc] initForReadingWithData:data];
        self.capturedView = [unarchiver decodeObjectForKey:kSavedView];
        [unarchiver finishDecoding];
        [unarchiver release];
        [data release];
    }
    else
        self.capturedView = [[[ModelViewCaptureContainer alloc] initWithSize:numMeshes/numFormats] autorelease];

    // Set to render the first mesh
    currRendererIndex = 0;
    meshFilename.text = [meshNames objectAtIndex:currRendererIndex];
    [self selectMeshAtIndex:currRendererIndex];

    if ([context API] != kEAGLRenderingAPIOpenGLES2)
    {
        // Create the lighting and material is using OpenGL ES 1.1
        GLfloat mat[4];
        mat[0]=0.2125;   mat[1]=0.1275;   mat[2]=0.054;   mat[3]=1.0;
        glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, mat);
        mat[0]=0.714;    mat[1]=0.4284;   mat[2]=0.18144;
        glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, mat);
        mat[0]=0.393548; mat[1]=0.271906; mat[2]=0.166721;
        GLfloat shine = 0.2;
        glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, mat);
        glMaterialf (GL_FRONT_AND_BACK, GL_SHININESS, shine * 128.0);

        GLfloat light_position [] = { 0.0,  0.0,   1.0, 0.0};
        glLightfv(GL_LIGHT0, GL_POSITION, light_position);

        glEnable(GL_LIGHTING);
        glEnable(GL_LIGHT0);
        glEnable(GL_NORMALIZE);
    }

	glEnable(GL_DEPTH_TEST);
	glEnable(GL_CULL_FACE);

	///// MobiMesh: end /////
}

- (void)dealloc
{
	if (program)
	{
		glDeleteProgram(program);
		program = 0;
	}

	// Tear down context.
	if ([EAGLContext currentContext] == context)
		[EAGLContext setCurrentContext:nil];
    
	[context release];

    if (vertexShaderSource != nil)
        [vertexShaderSource release];
    if (fragmentShaderSource != nil)
        [fragmentShaderSource release];

    for (int i=0; i<renderers.size(); ++i)
    {
        if (renderers[i] != nil)
            [renderers[i] release];
    }
    [capturedView release];
    [meshNames release];
    [meshFilename release];
    [errorMessageLabel release];
    [meshLoadingIndicator release];

	[super dealloc];

}

- (void)viewWillAppear:(BOOL)animated
{
	[self startAnimation];

	[super viewWillAppear:animated];
}

- (void)viewWillDisappear:(BOOL)animated
{
	[self stopAnimation];

	[super viewWillDisappear:animated];
}

- (void)viewDidUnload
{
	[super viewDidUnload];

	if (program)
	{
		glDeleteProgram(program);
		program = 0;
	}

	// Tear down context.
	if ([EAGLContext currentContext] == context)
		[EAGLContext setCurrentContext:nil];
	self.context = nil;	
}

- (NSInteger)animationFrameInterval
{
	return animationFrameInterval;
}

- (void)setAnimationFrameInterval:(NSInteger)frameInterval
{
	/*
	 Frame interval defines how many display frames must pass between each time the display link fires.
	 The display link will only fire 30 times a second when the frame internal is two on a display that refreshes 60 times a second. The default frame interval setting of one will fire 60 times a second when the display refreshes at 60 times a second. A frame interval setting of less than one results in undefined behavior.
	 */
	if (frameInterval >= 1)
	{
		animationFrameInterval = frameInterval;

		if (animating)
		{
			[self stopAnimation];
			[self startAnimation];
		}
	}
}

- (void)startAnimation
{
	if (!animating)
	{
		CADisplayLink *aDisplayLink = [CADisplayLink displayLinkWithTarget:self selector:@selector(drawFrame)];
		[aDisplayLink setFrameInterval:animationFrameInterval];
		[aDisplayLink addToRunLoop:[NSRunLoop currentRunLoop] forMode:NSDefaultRunLoopMode];
		self.displayLink = aDisplayLink;

		animating = TRUE;
	}
}

- (void)stopAnimation
{
	if (animating)
	{
		[self.displayLink invalidate];
		self.displayLink = nil;
		animating = FALSE;
	}
}

- (NSString*)dataFilePath
{
    NSArray* paths = NSSearchPathForDirectoriesInDomains(NSDocumentDirectory, NSUserDomainMask, YES);
    NSString* documentsDirectory = [paths objectAtIndex:0];
    return [documentsDirectory stringByAppendingPathComponent:kViewCaptureFilename];
}

- (void)saveCapturedData
{
    NSMutableData* data = [[NSMutableData alloc] init];
    NSKeyedArchiver* archiver = [[NSKeyedArchiver alloc] initForWritingWithMutableData:data];
    [archiver encodeObject:capturedView forKey:kSavedView];
    [archiver finishEncoding];
    [data writeToFile:[self dataFilePath] atomically:YES];
    [archiver release];
    [data release];
}

- (void)drawFrame
{
    // We have to check this flag before rendering because if the mesh is loading,
    // the OpenGL context does not belong to us and we cannot call OpenGL methods
    // in this case or OpenGL will crash.
    if (isMeshLoading)
        return;

	[(EAGLView *)self.view setFramebuffer];

    MobiMesh::IMeshRenderer* currRenderer = [self getCurrentRenderer];
    if (currRenderer != NULL)
    {
        if (autoRotate)
            currRenderer->rotate(1.0f, 0.0f, 1.0f, 0.0f);
        currRenderer->render();
    }
    else // invalid renderer, or mesh is still loading, just clear the screen
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	[(EAGLView *)self.view presentFramebuffer];
}

- (float)distanceBetweenTouches:(NSSet *)touches;
{
	int currentStage = 0;
	CGPoint point1, point2;
	
	for (UITouch *currentTouch in touches)
	{
		if (currentStage == 0)
		{
			point1 = [currentTouch locationInView:self.view];
			currentStage++;
		}
		else if (currentStage == 1) 
		{
			point2 = [currentTouch locationInView:self.view];
			currentStage++;
		}
	}
	return (sqrt((point1.x - point2.x) * (point1.x - point2.x) + (point1.y - point2.y) * (point1.y - point2.y)));
}

- (CGPoint)commonDirectionOfTouches:(NSSet *)touches;
{
	// Check to make sure that both fingers are moving in the same direction
	int currentStage = 0;
	CGPoint currentLocationOfTouch1, currentLocationOfTouch2, previousLocationOfTouch1, previousLocationOfTouch2;

	for (UITouch *currentTouch in touches)
	{
		if (currentStage == 0)
		{
			previousLocationOfTouch1 = [currentTouch previousLocationInView:self.view];
			currentLocationOfTouch1  = [currentTouch locationInView:self.view];
			currentStage++;
		}
		else if (currentStage == 1) 
		{
			previousLocationOfTouch2 = [currentTouch previousLocationInView:self.view];
			currentLocationOfTouch2  = [currentTouch locationInView:self.view];
			currentStage++;
		}
	}
	
	CGPoint directionOfTouch1, directionOfTouch2, commonDirection;
	// The sign of the Y touches is inverted, due to the inverted coordinate system of the iPhone
	directionOfTouch1.x =  currentLocationOfTouch1.x - previousLocationOfTouch1.x;
	directionOfTouch1.y = previousLocationOfTouch1.y -  currentLocationOfTouch1.y;
	directionOfTouch2.x =  currentLocationOfTouch2.x - previousLocationOfTouch2.x;
	directionOfTouch2.y = previousLocationOfTouch2.y -  currentLocationOfTouch2.y;

	// A two-finger movement should result in the direction of both touches being positive or negative
    // at the same time in X and Y
	if (!( ((directionOfTouch1.x <= 0) && (directionOfTouch2.x <= 0)) ||
          ((directionOfTouch1.x >= 0) && (directionOfTouch2.x >= 0)) ))
		return CGPointZero;
	if (!( ((directionOfTouch1.y <= 0) && (directionOfTouch2.y <= 0)) ||
          ((directionOfTouch1.y >= 0) && (directionOfTouch2.y >= 0)) ))
		return CGPointZero;

	// The movement ranges are averaged out
	commonDirection.x = ((directionOfTouch1.x + directionOfTouch1.x) / 2.0f);
	commonDirection.y = ((directionOfTouch1.y + directionOfTouch1.y) / 2.0f);
	
	return commonDirection;
}

- (void)touchesBegan:(NSSet *)touches withEvent:(UIEvent *)event
{
    NSMutableSet *currentTouches = [[[event touchesForView:self.view] mutableCopy] autorelease];
    [currentTouches minusSet:touches];
	
	// New touches are not yet included in the current touches for the view
	NSSet *totalTouches = [touches setByAddingObjectsFromSet:[event touchesForView:self.view]];
	if ([totalTouches count] > 1)
	{
		startingTouchDistance = [self distanceBetweenTouches:totalTouches];
		previousScale = 1.0f;
		twoFingersAreMoving  = NO;
		pinchGestureUnderway = NO;
		previousDirectionOfPanning = CGPointZero;
	}
	else
	{
		lastMovementPosition = [[touches anyObject] locationInView:self.view];
	}
    // *** MobiMesh Tutorial: Double tap to reset camera
    UITouch* touch = [touches anyObject];
    NSUInteger numTaps = [touch tapCount];
    RendererState* currRendererState = renderers[currRendererIndex];
    if (numTaps > 1 && currRendererState != nil)
    {
        MobiMesh::IMeshRenderer* currRenderer = rendererFactory.get_renderer(currRendererState.rendererHandle);
        if (currRenderer)
            currRenderer->set_default_view();
    }
}

- (void)touchesMoved:(NSSet *)touches withEvent:(UIEvent *)event
{	
	if ([[event touchesForView:self.view] count] > 1) // Pinch gesture, possibly two-finger movement
	{
		CGPoint directionOfPanning = CGPointZero;
		
		// Two finger panning
		if ([touches count] > 1) // Check to make sure that both fingers are moving
		{
			directionOfPanning = [self commonDirectionOfTouches:touches];
		}
		
		if ( (directionOfPanning.x != 0) || (directionOfPanning.y != 0) ) // Don't scale while doing the two-finger panning
		{
			if (pinchGestureUnderway)
			{
				if (sqrt(previousDirectionOfPanning.x * previousDirectionOfPanning.x + previousDirectionOfPanning.y * previousDirectionOfPanning.y) > 0.1 )
				{
					pinchGestureUnderway = NO;
				}
				previousDirectionOfPanning.x += directionOfPanning.x;
				previousDirectionOfPanning.y += directionOfPanning.y;
			}
			if (!pinchGestureUnderway)
			{
				twoFingersAreMoving = YES;
                // Translate by an amount relative to the size of the mesh
                RendererState* currRendererState = renderers[currRendererIndex];
                if (currRendererState != nil)
                {
                    MobiMesh::IMeshRenderer* renderer = rendererFactory.get_renderer(currRendererState.rendererHandle);
                    if (renderer)
                    {
                        float xPanning = directionOfPanning.x / ((EAGLView *)self.view).framebufferWidth *
                            renderer->get_bounding_length();
                        float yPanning = directionOfPanning.y / ((EAGLView *)self.view).framebufferHeight *
                            renderer->get_bounding_length();
                        renderer->translate(xPanning, yPanning, 0.0);
                    }
                }
				[self drawFrame];
				previousDirectionOfPanning = CGPointZero;
			}
		}
		else
		{
			float newTouchDistance = [self distanceBetweenTouches:[event touchesForView:self.view]];
			if (twoFingersAreMoving)
			{
				// If fingers have moved more than 10% apart, start pinch gesture again
				if ( fabs(1 - (newTouchDistance / startingTouchDistance) / previousScale) > 0.3 )
				{
					twoFingersAreMoving = NO;
				}
			}
			if (!twoFingersAreMoving)
			{
				// Scale using pinch gesture
                RendererState* currRendererState = renderers[currRendererIndex];
                if (currRendererState != nil)
                {
                    MobiMesh::IMeshRenderer* currRenderer = rendererFactory.get_renderer(currRendererState.rendererHandle);
                    if (currRenderer)
                        currRenderer-> scale((newTouchDistance / startingTouchDistance) / previousScale);
                }
				[self drawFrame];
				previousScale = (newTouchDistance / startingTouchDistance);
				pinchGestureUnderway = YES;
			}
		}
	}	
	else // Single-touch rotation of object
	{
		CGPoint currentMovementPosition = [[touches anyObject] locationInView:self.view];

        RendererState* currRendererState = renderers[currRendererIndex];
        if (currRendererState != nil)
        {
            MobiMesh::IMeshRenderer* renderer = rendererFactory.get_renderer(currRendererState.rendererHandle);
            if (renderer)
            {
                MobiMesh::TrackBall trackBall;
                float angle_rad, axis_x, axis_y, axis_z;
                bool rotate = trackBall.get_rotation(lastMovementPosition.x,
                                                     ((EAGLView *)self.view).framebufferHeight - lastMovementPosition.y,
                                                     currentMovementPosition.x,
                                                     ((EAGLView *)self.view).framebufferHeight - currentMovementPosition.y,
                                                     ((EAGLView *)self.view).framebufferWidth,
                                                     ((EAGLView *)self.view).framebufferHeight,
                                                     angle_rad, axis_x, axis_y, axis_z);
                if (rotate)
                    renderer->rotate(angle_rad * 180.0 / M_PI, axis_x, axis_y, axis_z);
            }
        }

		[self drawFrame];
		lastMovementPosition = currentMovementPosition;
	}
}

- (void)touchesCancelled:(NSSet *)touches withEvent:(UIEvent *)event 
{
	NSMutableSet *remainingTouches = [[[event touchesForView:self.view] mutableCopy] autorelease];
    [remainingTouches minusSet:touches];
	if ([remainingTouches count] < 2)
	{
		twoFingersAreMoving  = NO;
		pinchGestureUnderway = NO;
		previousDirectionOfPanning = CGPointZero;
		
		lastMovementPosition = [[remainingTouches anyObject] locationInView:self.view];
	}
}

- (MobiMesh::IMeshRenderer*)getCurrentRenderer
{
    return [self getRendererAtIndex:currRendererIndex];
}

- (MobiMesh::IMeshRenderer*)getRendererAtIndex:(int)rendererIndex
{
    RendererState* currRendererState = renderers[rendererIndex];
    if (!currRendererState)
        return NULL;
    // Return pointer to renderer, or NULL if handle is invalid
    return rendererFactory.get_renderer(currRendererState.rendererHandle);
}

- (IBAction)showMenu:(id)sender
{
    // If mesh is loading, disable UI
    if (isMeshLoading)
        return;

    MobiMeshAppDelegate* appDelegate = (MobiMeshAppDelegate*)[[UIApplication sharedApplication] delegate];
    [appDelegate showMenu];
}

- (IBAction)toggleAutoRotate:(id)sender
{
    if (isMeshLoading)
        return;

    autoRotate = !autoRotate;
}

- (IBAction)captureView:(id)sender
{
    // If mesh is loading, disable UI
    if (isMeshLoading)
        return;

    // Do nothing if renderer is invalid
    if ([self getCurrentRenderer] == NULL)
        return;

    [self captureCurrentView];

    UIAlertView* alert = [[UIAlertView alloc]
                          initWithTitle:@"Capture"
                          message:@"View captured"
                          delegate:nil
                          cancelButtonTitle:@"OK"
                          otherButtonTitles:nil, nil];
    [alert show];
    [alert release];
}

- (IBAction)loadView:(id)sender
{
    // If mesh is loading, disable UI
    if (isMeshLoading)
        return;

    MobiMesh::IMeshRenderer* currRenderer = [self getCurrentRenderer];
    if (!currRenderer)
        return;

    MobiMeshAppDelegate* appDelegate = (MobiMeshAppDelegate*)[[UIApplication sharedApplication] delegate];
    [appDelegate showCapturesMenu:currRendererIndex/numFormats];
}

- (IBAction)previousMesh:(id)sender
{
    // If mesh is loading, disable UI
    if (isMeshLoading)
        return;

    int rendererIndex = currRendererIndex;
    if (rendererIndex == 0)
    {
        // Loop back to last mesh
        rendererIndex = [meshNames count];
    }
    --rendererIndex;

    [self selectMeshAtIndex:rendererIndex];
}

- (IBAction)nextMesh:(id)sender
{
    // If mesh is loading, disable UI
    if (isMeshLoading)
        return;

    int rendererIndex = currRendererIndex+1;
    if (rendererIndex == [meshNames count])
        rendererIndex = 0;

    [self selectMeshAtIndex:rendererIndex];
}

- (void)selectMeshAtIndex:(int)rendererIndex
{
    if (rendererIndex >= [meshNames count] || rendererIndex < 0)
        return;

    int prevRendererIndex = currRendererIndex;
    currRendererIndex = rendererIndex;

    [meshFilename setText:[meshNames objectAtIndex:currRendererIndex]];
    [errorMessageLabel setText:nil]; // clear any error messages from screen

    int prevMeshIndex = prevRendererIndex/numFormats;
    int currMeshIndex = currRendererIndex/numFormats;
    if (prevMeshIndex != currMeshIndex)
    {
        // If we are displaying a new mesh (i.e. not just different format,
        // but different mesh altogether), unload the previous mesh[es]
        [self unloadMeshAtIndex:(3*prevMeshIndex)];
        [self unloadMeshAtIndex:(3*prevMeshIndex+1)];
        [self unloadMeshAtIndex:(3*prevMeshIndex+2)];
    }
    else if (prevRendererIndex != currRendererIndex)
    {
        // Otherwise, deactivate the previous mesh to save some memory,
        // but do not unload it. This is to allow users to cycle between
        // the formats quickly without having to reload and re-apply the
        // refinements.
        MobiMesh::IMeshRenderer* prevRenderer =
            [self getRendererAtIndex:prevRendererIndex];
        if (prevRenderer)
            prevRenderer->deactivate();
    }

    [self loadCurrMeshIfNotLoaded];
}

- (void)unloadMeshAtIndex:(int)rendererIndex
{
    RendererState* currRendererState = renderers[rendererIndex];
    if (currRendererState)
        [currRendererState unloadRenderer];
}

- (void)loadMeshAtIndex:(int)index
{
    if (index >= [meshNames count]) // invalid index
        return;

    NSString* appFolderPath = [[[NSBundle mainBundle] resourcePath] stringByAppendingString:@"/"];
    NSString* meshPath = [appFolderPath stringByAppendingString:[meshNames objectAtIndex:index]];

    RendererState* currRendererState = renderers[index];
    if (!currRendererState)
    {
        // create the renderer state for the first time and load the mesh.
        if ([context API] == kEAGLRenderingAPIOpenGLES2)
        {
            currRendererState = [[RendererState alloc] initWithFactory: &rendererFactory
                                                          withFilename: meshPath
                                                      withVertexShader: vertexShaderSource
                                                    withFragmentShader: fragmentShaderSource];
        }
        else
        {
            currRendererState = [[RendererState alloc] initWithFactory: &rendererFactory
                                                          withFilename: meshPath];
        }

        if (currRendererState == nil)
        {
            // Failed to alloc memory for renderer state
            [self performSelectorOnMainThread:@selector(setErrorString:)
                                   withObject:@"Unable to create renderer state"
                                waitUntilDone:NO];
            return;
        }

        renderers[index] = currRendererState;
    }

    // Attempt to load the renderer/mesh
    try
    {
        [currRendererState loadRendererWithViewportWidth:((EAGLView *)self.view).framebufferWidth
                                      withViewportHeight:((EAGLView *)self.view).framebufferHeight];
        MobiMesh::IMeshRenderer* currRenderer = [self getRendererAtIndex:index];

        if (currRenderer == NULL)
        {
            // If current renderer is NULL despite having loaded the renderer, then the
            // factory is unable to create more renderers, indicate this error
            [self performSelectorOnMainThread:@selector(setErrorString:)
                                   withObject:@"Too many active renderers, "
                                               "unload some first."
                                waitUntilDone:NO];
        }
    }
    catch (MobiMesh::InvalidMeshTypeException& ex)
    {
        [self performSelectorOnMainThread:@selector(setErrorString:)
                               withObject:[NSString stringWithFormat:@"%s", ex.what()]
                            waitUntilDone:NO];
    }
    catch (MobiMesh::UnableToCreateRendererException& ex)
    {
        [self performSelectorOnMainThread:@selector(setErrorString:)
                               withObject:[NSString stringWithFormat:@"%s", ex.what()]
                            waitUntilDone:NO];
    }
}

- (void) setErrorString:(NSString *)error
{
    [errorMessageLabel setText:error];
}

- (void)captureCurrentView
{
    [self captureViewAtRendererIndex:currRendererIndex];
}

- (void)captureViewAtRendererIndex:(int)rendererIndex
{
    MobiMesh::IMeshRenderer* currRenderer = [self getRendererAtIndex:rendererIndex];
    if (!currRenderer)
        return;

    // Capture the default view when the mesh is loaded
    std::pair<glm::mat4, glm::mat4> matrix;	
    currRenderer->capture_view(matrix.first, matrix.second);
    [capturedView addModelViewMatrix:matrix atMeshIndex:rendererIndex/numFormats];
}

- (void)loadViewAtViewIndex:(int)viewIndex
{
    MobiMesh::IMeshRenderer* currRenderer = [self getCurrentRenderer];
    if (!currRenderer)
        return;

    std::pair<glm::mat4, glm::mat4> matrix =
        [capturedView getModelViewMatrix:currRendererIndex/numFormats atViewIndex:viewIndex];
    currRenderer->set_view(matrix.first, matrix.second);
}

// Thread for loading mesh in background
- (void)loadMeshInBackground
{
    NSAutoreleasePool* pool = [[NSAutoreleasePool alloc] init];
    
    // Set the current context of this thread to the context we have stored
    // This is necessary because loading the mesh/renderers make calls to OpenGL
    // and the context must be set in order for these calls to work correctly.
    [EAGLContext setCurrentContext:context];

    // Load the renderer if necessary
    RendererState* currRendererState = renderers[currRendererIndex];
    if (!currRendererState || ![currRendererState isLoaded])
        [self loadMeshAtIndex:currRendererIndex];

    // Activate the renderer if it is inactive
    MobiMesh::IMeshRenderer* currRenderer = [self getCurrentRenderer];
    if (currRenderer && !currRenderer->is_active())
        currRenderer->activate();

    [self performSelectorOnMainThread:@selector(didStopMeshLoading) withObject:nil waitUntilDone:NO];

    [pool release];
}

- (void)loadCurrMeshIfNotLoaded
{
    RendererState* currRendererState = renderers[currRendererIndex];
    MobiMesh::IMeshRenderer* currRenderer = [self getCurrentRenderer];

    if (!currRendererState || ![currRendererState isLoaded] ||
        !currRenderer->is_active())
    {
        // Renderer has not been loaded or is inactive,
        // load it in the background first
        [self willStartMeshLoading];
        [self performSelectorInBackground:@selector(loadMeshInBackground) withObject:nil];
    }
}

- (void)willStartMeshLoading
{
    [meshLoadingIndicator startAnimating];

    // Clear the screen in case there is currently another mesh on it.
	[(EAGLView *)self.view setFramebuffer];    
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	[(EAGLView *)self.view presentFramebuffer];

    // Set this flag to disable the UI so that we do not load multiple
    // meshes at the same time
    isMeshLoading = true;
}

- (void)didStopMeshLoading
{
    [meshLoadingIndicator stopAnimating];

    // Set the OpenGL context for this thread (since loading the mesh
    // possibly transfered this context to the mesh loading thread)
    [(EAGLView *)self.view setContext:context];

    // Clear this flag so that UI can work again
    isMeshLoading = false;

    // Redraw before we return because this method can take a while if it
    // needs to refine a view dependent mesh. This method upon returning will
    // stop spinning the meshLoadingIndicator. If we do not call drawFrame
    // here, users will see a black screen without the indicator until the
    // next frame is drawn, which is not desirable if it takes very long to
    // draw a frame.
    [self drawFrame];
}

///// MobiMesh: end /////

- (void)didReceiveMemoryWarning
{
	// Releases the view if it doesn't have a superview.
	[super didReceiveMemoryWarning];

	// Release any cached data, images, etc. that aren't in use.
    // Unload all meshes except the current one
    for (int i=0; i<[meshNames count]; ++i)
        if (i != currRendererIndex)
            [self unloadMeshAtIndex:i];
}

@end
