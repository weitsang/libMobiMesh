//
//  TVPMRendererAppDelegate.h
//  TVPMRenderer
//
//  Created by CS4344 on 6/30/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//

#import <UIKit/UIKit.h>

@class TVPMRendererViewController;

@interface TVPMRendererAppDelegate : NSObject <UIApplicationDelegate> {

}

@property (nonatomic, retain) IBOutlet UIWindow *window;

@property (nonatomic, retain) IBOutlet TVPMRendererViewController *viewController;

@end
