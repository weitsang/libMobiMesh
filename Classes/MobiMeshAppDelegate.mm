//
//  MobiMeshAppDelegate.m
//  MobiMesh
//
//  Created by Ong Yuh Shin on 3/4/11.
//  Copyright 2011 National University of Singapore. All rights reserved.
//

#import "MobiMeshAppDelegate.h"
#import "MobiMeshViewController.h"
#import "MenuViewController.h"
#import "ModelviewManagerViewController.h"

@implementation MobiMeshAppDelegate

@synthesize window;
@synthesize meshViewController, menuViewController, modelviewManagerViewController;

- (BOOL)application:(UIApplication *)application didFinishLaunchingWithOptions:(NSDictionary *)launchOptions
{
    menuViewController.meshNames = meshViewController.meshNames;
    modelviewManagerViewController.captures = meshViewController.capturedView;

    [self.window addSubview:self.meshViewController.view];
    return YES;
}

- (void)applicationWillResignActive:(UIApplication *)application
{
    [self.meshViewController stopAnimation];
}

- (void)applicationDidBecomeActive:(UIApplication *)application
{
    [self.meshViewController startAnimation];
}

- (void)applicationWillTerminate:(UIApplication *)application
{
    [self.meshViewController stopAnimation];
    [self.meshViewController saveCapturedData];
}

- (void)applicationDidEnterBackground:(UIApplication *)application
{
    // Handle any background procedures not related to animation here.
}

- (void)applicationWillEnterForeground:(UIApplication *)application
{
    // Handle any foreground procedures not related to animation here.
}

- (void)dealloc
{
    [window release];
    [meshViewController release];
    [menuViewController release];
    [modelviewManagerViewController release];

    [super dealloc];
}

- (IBAction)showMenu
{
    [menuViewController.meshTable reloadData];
    [self.window addSubview:self.menuViewController.view];
    [meshViewController.view removeFromSuperview];
}

- (void)showMesh:(int)rendererIndex
{
    [self showMesh];
    [meshViewController selectMeshAtIndex:rendererIndex];
}

- (void)showMeshAtMeshIndex:(int)meshIndex atViewIndex:(int)viewIndex
{
    [self showMesh];
    [meshViewController loadViewAtViewIndex:viewIndex];
}

- (void)showMesh
{
    [menuViewController.view removeFromSuperview];
    [modelviewManagerViewController.view removeFromSuperview];
    [self.window addSubview:self.meshViewController.view];
}

- (void)showCapturesMenu:(int)meshIndex
{
    modelviewManagerViewController.meshIndex = meshIndex;
    [modelviewManagerViewController.captureTable reloadData];

    [self.window addSubview:self.modelviewManagerViewController.view];
    [meshViewController.view removeFromSuperview];
}

@end
