//
//  ViewManagerViewController.h
//  MobiMesh
//
//  Created by Ong Yuh Shin on 5/23/11.
//  Copyright 2011 National University of Singapore. All rights reserved.
//

#import <UIKit/UIKit.h>

#import "ModelViewCaptureContainer.h"


@interface ModelviewManagerViewController : UIViewController <UITableViewDelegate, UITableViewDataSource> {
    ModelViewCaptureContainer* captures;
    IBOutlet UITableView* captureTable;
    IBOutlet UIBarButtonItem* deleteModeToggleButton;
    int meshIndex; // to indicate which mesh's captures to display
}

@property (retain, nonatomic) ModelViewCaptureContainer* captures;
@property (retain, nonatomic) UITableView* captureTable;
@property (retain, nonatomic) UIBarButtonItem* deleteModeToggleButton;
@property (nonatomic) int meshIndex;

- (NSInteger)tableView:(UITableView *)tableView numberOfRowsInSection:(NSInteger)section;
- (UITableViewCell*)tableView:(UITableView *)tableView cellForRowAtIndexPath:(NSIndexPath *)indexPath;
- (void)tableView:(UITableView *)tableView didSelectRowAtIndexPath:(NSIndexPath *)indexPath;
- (UITableViewCellEditingStyle)tableView:(UITableView *)tableView editingStyleForRowAtIndexPath:(NSIndexPath *)indexPath;
- (void)tableView:(UITableView *)tableView commitEditingStyle:(UITableViewCellEditingStyle)editingStyle forRowAtIndexPath:(NSIndexPath *)indexPath;

- (IBAction)back:(id)sender;

@end
