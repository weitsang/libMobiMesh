//
//  MenuViewController.h
//  MobiMesh
//
//  Created by Ong Yuh Shin on 5/11/11.
//  Copyright 2011 National University of Singapore. All rights reserved.
//

#import <UIKit/UIKit.h>
#import <vector>


@interface MenuViewController : UIViewController <UITableViewDelegate, UITableViewDataSource> {
    NSArray* meshNames;
    IBOutlet UITableView* meshTable;

    int currRendererIndex;
}

@property (retain, nonatomic) NSArray* meshNames;
@property (retain, nonatomic) UITableView* meshTable;

- (NSInteger)tableView:(UITableView *)tableView numberOfRowsInSection:(NSInteger)section;
- (UITableViewCell*)tableView:(UITableView *)tableView cellForRowAtIndexPath:(NSIndexPath *)indexPath;
- (void)tableView:(UITableView *)tableView didSelectRowAtIndexPath:(NSIndexPath *)indexPath;

- (IBAction)back:(id)sender;

@end
