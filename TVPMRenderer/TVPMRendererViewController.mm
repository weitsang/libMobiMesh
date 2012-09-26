//
//  TVPMRendererViewController.m
//  MobiMesh
//
//  Created by Ong Yuh Shin on 6/30/11.
//  Copyright 2011 National University of Singapore. All rights reserved.
//

#import <QuartzCore/QuartzCore.h>

#import "TVPMRendererViewController.h"
#import "EAGLView.h"

#import "GLMeshRendererES1.h"
#import "GLMeshRendererES2.h"
#import "TrulyViewDependentPM32.h"
#import "FileReadErrorException.h"
#import "InvalidFileFormatException.h"
#import "ShaderCompileException.h"
#import "ProgramLinkingException.h"
#import "InvalidGLProgramException.h"

// Uniform index.
enum {
    UNIFORM_TRANSLATE,
    NUM_UNIFORMS
};
GLint uniforms[NUM_UNIFORMS];

// Attribute index.
enum {
    ATTRIB_VERTEX,
    ATTRIB_COLOR,
    NUM_ATTRIBUTES
};

@interface TVPMRendererViewController ()
@property (nonatomic, retain) EAGLContext *context;
@property (nonatomic, assign) CADisplayLink *displayLink;
@end

@implementation TVPMRendererViewController

@synthesize animating, context, displayLink;

- (void)awakeFromNib
{
    EAGLContext *aContext = [[EAGLContext alloc] initWithAPI:kEAGLRenderingAPIOpenGLES2];
    
    if (!aContext) {
        aContext = [[EAGLContext alloc] initWithAPI:kEAGLRenderingAPIOpenGLES1];
    }
    
    if (!aContext)
        NSLog(@"Failed to create ES context");
    else if (![EAGLContext setCurrentContext:aContext])
        NSLog(@"Failed to set ES context current");
    
	self.context = aContext;
	[aContext release];
	
    [(EAGLView *)self.view setContext:context];
    [(EAGLView *)self.view setFramebuffer];

    // *** MobiMesh Tutorial: Create the mesh
    NSString *meshPath = [[NSBundle mainBundle] pathForResource:@"thai" ofType:@"tvpm"];
    try
    {
        mesh = new MobiMesh::TrulyViewDependentPM32([meshPath cStringUsingEncoding:NSASCIIStringEncoding]);
    }
    catch (MobiMesh::FileReadErrorException& ex)
    {
        NSLog(@"%s", ex.what());
        exit(1);
    }
    catch (MobiMesh::InvalidFileFormatException& ex)
    {
        NSLog(@"%s", ex.what());
        exit(1);
    }

    if ([context API] == kEAGLRenderingAPIOpenGLES2)
    {
        // *** MobiMesh Tutorial: Load shaders for renderer if using ES2
        NSString *vertexShaderFile = [[NSBundle mainBundle] pathForResource:@"Shader" ofType:@"vsh"];
        NSString *fragmentShaderFile = [[NSBundle mainBundle] pathForResource:@"Shader" ofType:@"fsh"];
        NSString* vertexShaderSource =
            [NSString stringWithContentsOfFile:vertexShaderFile
                                      encoding:NSUTF8StringEncoding error:nil];        
        NSString* fragmentShaderSource =
            [NSString stringWithContentsOfFile:fragmentShaderFile
                                      encoding:NSUTF8StringEncoding error:nil];
        if (!fragmentShaderSource || !vertexShaderSource)
        {
            NSLog(@"Failed to load shaders");
            exit(1);
        }

        // *** MobiMesh Tutorial: Create the renderer and exit on failure.
        try
        {
            renderer = new MobiMesh::GLMeshRendererES2([vertexShaderSource cStringUsingEncoding:NSUTF8StringEncoding],
                                                       [fragmentShaderSource cStringUsingEncoding:NSUTF8StringEncoding],
                                                       [(EAGLView *)self.view framebufferWidth],
                                                       [(EAGLView *)self.view framebufferHeight]);
        }
        catch (MobiMesh::ShaderCompileException& ex)
        {
            NSLog(@"%s", ex.what());
            exit(1);
        }
        catch (MobiMesh::ProgramLinkingException& ex)
        {
            NSLog(@"%s", ex.what());
            exit(1);
        }
        catch (MobiMesh::InvalidGLProgramException& ex)
        {
            NSLog(@"%s", ex.what());
            exit(1);
        }
    }
    else
    {
        // *** MobiMesh Tutorial: Create ES1 renderer if using ES1
        renderer = new MobiMesh::GLMeshRendererES1([(EAGLView *)self.view framebufferWidth],
                                                   [(EAGLView *)self.view framebufferHeight]);
        
        // *** MobiMesh Tutorial: Set up lighting and material
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

    // *** MobiMesh Tutorial: Create the TVPM renderer
    tvpmRenderer = new MobiMesh::GLTvpmMeshRenderer(*renderer, mesh);

    // *** MobiMesh Tutorial: Enable depth testing for more accurate rendering
	glEnable(GL_DEPTH_TEST);
    // *** MobiMesh Tutorial: Enable back face culling
	glEnable(GL_CULL_FACE);

    animating = FALSE;
    animationFrameInterval = 1;
    self.displayLink = nil;
}

- (void)dealloc
{
    if (program) {
        glDeleteProgram(program);
        program = 0;
    }
    
    // Tear down context.
    if ([EAGLContext currentContext] == context)
        [EAGLContext setCurrentContext:nil];
    
    [context release];

    // *** MobiMesh Tutorial: Release renderers and mesh
    if (tvpmRenderer)
        delete tvpmRenderer;
    if (renderer)
        delete renderer;
    if (mesh)
        delete mesh;

    [super dealloc];
}

- (void)didReceiveMemoryWarning
{
    // Releases the view if it doesn't have a superview.
    [super didReceiveMemoryWarning];
    
    // Release any cached data, images, etc. that aren't in use.
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
	
    if (program) {
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
    if (frameInterval >= 1) {
        animationFrameInterval = frameInterval;
        
        if (animating) {
            [self stopAnimation];
            [self startAnimation];
        }
    }
}

- (void)startAnimation
{
    if (!animating) {
        CADisplayLink *aDisplayLink = [[UIScreen mainScreen] displayLinkWithTarget:self selector:@selector(drawFrame)];
        [aDisplayLink setFrameInterval:animationFrameInterval];
        [aDisplayLink addToRunLoop:[NSRunLoop currentRunLoop] forMode:NSDefaultRunLoopMode];
        self.displayLink = aDisplayLink;
        
        animating = TRUE;
    }
}

- (void)stopAnimation
{
    if (animating) {
        [self.displayLink invalidate];
        self.displayLink = nil;
        animating = FALSE;
    }
}

- (void)drawFrame
{
    [(EAGLView *)self.view setFramebuffer];

    // *** MobiMesh Tutorial: Rotate by 1 degree about y-axis in each frame and render the mesh
    tvpmRenderer->rotate(1.0f, 0.0f, 1.0f, 0.0f);
    tvpmRenderer->render();

    [(EAGLView *)self.view presentFramebuffer];
}

@end
