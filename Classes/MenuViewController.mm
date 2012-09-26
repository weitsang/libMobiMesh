//
//  MenuViewController.m
//  MobiMesh
//
//  Created by Ong Yuh Shin on 5/11/11.
//  Copyright 2011 National University of Singapore. All rights reserved.
//

#import "MenuViewController.h"
#import "MobiMeshAppDelegate.h"


@implementation MenuViewController

@synthesize meshTable;
@synthesize meshNames;

- (id)initWithNibName:(NSString *)nibNameOrNil bundle:(NSBundle *)nibBundleOrNil
{
    self = [super initWithNibName:nibNameOrNil bundle:nibBundleOrNil];
    if (self) {
        // Custom initialization
    }
    return self;
}

- (void)dealloc
{
    [meshTable release];
    [meshNames release];
    [super dealloc];
}

- (void)didReceiveMemoryWarning
{
    // Releases the view if it doesn't have a superview.
    [super didReceiveMemoryWarning];
    
    // Release any cached data, images, etc that aren't in use.
}

#pragma mark - View lifecycle

- (void)viewDidLoad
{
    [super viewDidLoad];
    // Do any additional setup after loading the view from its nib.

    meshTable.backgroundColor = [UIColor grayColor];
    meshTable.dataSource = self;
    meshTable.delegate = self;
}

- (void)viewDidUnload
{
    [super viewDidUnload];
    // Release any retained subviews of the main view.
    // e.g. self.myOutlet = nil;
}

- (BOOL)shouldAutorotateToInterfaceOrientation:(UIInterfaceOrientation)interfaceOrientation
{
    // Return YES for supported orientations
    return (interfaceOrientation == UIInterfaceOrientationPortrait);
}

- (NSInteger)tableView:(UITableView *)tableView numberOfRowsInSection:(NSInteger)section
{
    return [meshNames count];
}

- (UITableViewCell*)tableView:(UITableView *)tableView cellForRowAtIndexPath:(NSIndexPath *)indexPath
{
    static NSString* identifier = @"MobiMeshId";
    
    UITableViewCell* cell = [tableView dequeueReusableCellWithIdentifier:identifier];
    if (cell == nil)
    {
        cell = [[[UITableViewCell alloc] initWithFrame:CGRectZero reuseIdentifier:identifier] autorelease];
        cell.accessoryType = UITableViewCellAccessoryDisclosureIndicator;
        cell.textLabel.textColor = [UIColor whiteColor];
        cell.textLabel.shadowColor = [UIColor blackColor];
        cell.textLabel.shadowOffset = CGSizeMake(1, 1);
        cell.selectionStyle = UITableViewCellSelectionStyleGray;
    }
    
    NSString* meshName = [meshNames objectAtIndex:indexPath.row];
    [cell.textLabel setText:meshName];

    return cell;
}

- (void)tableView:(UITableView *)tableView didSelectRowAtIndexPath:(NSIndexPath *)indexPath
{
    MobiMeshAppDelegate* appDelegate = (MobiMeshAppDelegate*)[[UIApplication sharedApplication] delegate];
    [appDelegate showMesh:indexPath.row];
}

- (IBAction)back:(id)sender
{
    MobiMeshAppDelegate* appDelegate = (MobiMeshAppDelegate*)[[UIApplication sharedApplication] delegate];
    [appDelegate showMesh];
}

@end
