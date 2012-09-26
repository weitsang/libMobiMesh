//
//  MobiMeshAppDelegate.h
//  MobiMesh
//
//  Created by Ong Yuh Shin on 3/4/11.
//  Copyright 2011 National University of Singapore. All rights reserved.
//

#import <UIKit/UIKit.h>

@class MobiMeshViewController;
@class MenuViewController;
@class ModelviewManagerViewController;

@interface MobiMeshAppDelegate : NSObject <UIApplicationDelegate> {
    UIWindow *window;
    MobiMeshViewController *meshViewController;
    MenuViewController *menuViewController;
    ModelviewManagerViewController *modelviewManagerViewController;
}

@property (nonatomic, retain) IBOutlet UIWindow *window;
@property (nonatomic, retain) IBOutlet MobiMeshViewController *meshViewController;
@property (nonatomic, retain) IBOutlet MenuViewController *menuViewController;
@property (nonatomic, retain) IBOutlet ModelviewManagerViewController *modelviewManagerViewController;

- (void)showMenu;
- (void)showMesh:(int)rendererIndex;
- (void)showMeshAtMeshIndex:(int)meshIndex atViewIndex:(int)viewIndex;
- (void)showMesh;
- (void)showCapturesMenu:(int)meshIndex;

@end

