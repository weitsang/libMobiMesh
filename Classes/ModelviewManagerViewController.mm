//
//  ViewManagerViewController.m
//  MobiMesh
//
//  Created by Ong Yuh Shin on 5/23/11.
//  Copyright 2011 National University of Singapore. All rights reserved.
//

#import "MobiMeshAppDelegate.h"
#import "ModelviewManagerViewController.h"


@implementation ModelviewManagerViewController

@synthesize captures, captureTable, meshIndex, deleteModeToggleButton;

- (id)initWithNibName:(NSString *)nibNameOrNil bundle:(NSBundle *)nibBundleOrNil
{
    self = [super initWithNibName:nibNameOrNil bundle:nibBundleOrNil];
    if (self) {
        // Custom initialization
        self.meshIndex = 0;
    }
    return self;
}

- (void)dealloc
{
    [captures release];
    [captureTable release];
    [deleteModeToggleButton release];

    [super dealloc];
}

- (void)didReceiveMemoryWarning
{
    // Releases the view if it doesn't have a superview.
    [super didReceiveMemoryWarning];
    
    // Release any cached data, images, etc that aren't in use.
}

#pragma mark - View lifecycle

/*
// Implement loadView to create a view hierarchy programmatically, without using a nib.
- (void)loadView
{
}
*/

// Implement viewDidLoad to do additional setup after loading the view, typically from a nib.
- (void)viewDidLoad
{
    [super viewDidLoad];

    captureTable.backgroundColor = [UIColor grayColor];
    captureTable.dataSource = self;
    captureTable.delegate = self;
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
    return [captures numViews:meshIndex];
}

- (UITableViewCell*)tableView:(UITableView *)tableView cellForRowAtIndexPath:(NSIndexPath *)indexPath
{
    static NSString* identifier = @"MobiMeshCapturesId";

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

    NSString* cellText = [[NSString alloc]initWithFormat:@"Capture %d", indexPath.row];
    [cell.textLabel setText:cellText];
    [cellText release];

    return cell;
}

- (void)tableView:(UITableView *)tableView didSelectRowAtIndexPath:(NSIndexPath *)indexPath
{
    if (!captureTable.editing)
    {
        MobiMeshAppDelegate* appDelegate = (MobiMeshAppDelegate*)[[UIApplication sharedApplication] delegate];
        [appDelegate showMeshAtMeshIndex:meshIndex atViewIndex:indexPath.row];
    }
}

- (UITableViewCellEditingStyle)tableView:(UITableView *)tableView editingStyleForRowAtIndexPath:(NSIndexPath *)indexPath
{
    if (captureTable.editing)
        return UITableViewCellEditingStyleDelete;
    else
        return UITableViewCellEditingStyleNone;
}

- (void)tableView:(UITableView *)tableView commitEditingStyle:(UITableViewCellEditingStyle)editingStyle forRowAtIndexPath:(NSIndexPath *)indexPath
{
    if (editingStyle != UITableViewCellEditingStyleDelete)
        return;

    [captures removeModelViewMatrixAtMeshIndex:meshIndex atViewIndex:indexPath.row];
    [captureTable reloadData];
}

- (IBAction)back:(id)sender
{
    MobiMeshAppDelegate* appDelegate = (MobiMeshAppDelegate*)[[UIApplication sharedApplication] delegate];
    [appDelegate showMesh];
}

- (IBAction)toggleDeleteMode:(id)sender
{
    [captureTable setEditing:!captureTable.editing animated:YES];

    if (captureTable.editing)
    {
        NSString* btnTitle = [[NSString alloc] initWithFormat:@"Done"];
        [deleteModeToggleButton setTitle:btnTitle];
        [btnTitle release];
    }
    else
    {
        NSString* btnTitle = [[NSString alloc] initWithFormat:@"Delete"];
        [deleteModeToggleButton setTitle:btnTitle];
        [btnTitle release];
    }
}

@end
