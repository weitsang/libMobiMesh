//
//  ModelViewCaptureContainer.h
//  MobiMesh
//
//  Created by Ong Yuh Shin on 5/23/11.
//  Copyright 2011 National University of Singapore. All rights reserved.
//

#import <Foundation/Foundation.h>

#import <vector>

#import "glm.hpp"


@interface ModelViewCaptureContainer : NSObject<NSCoding, NSCopying> {
    // A doubly indexed/2D vector. A modelview matrix is made up of the pair
    // (lookat matrix, transformation matrix). The actual matrix can be
    // obtained from lookat * transformation. To access the actual matrix,
    // we use the meshIndex and viewIndex:
    // modelViewMatrices[meshIndex][viewIndex]
    // The meshIndex indicates which mesh the matrix describes, and the
    // viewIndex
    // indicates the particular view for the particular mesh identified by
    // the meshIndex.
    std::vector<std::vector<std::pair<glm::mat4, glm::mat4> > > modelViewMatrices;
}

- (id)initWithSize:(int)size;
- (id)initWithCoder:(NSCoder *)decoder;

- (BOOL)setModelViewMatrix:(std::pair<glm::mat4, glm::mat4>)matrix atMeshIndex:(int)meshIndex atViewIndex:(int)viewIndex;
- (int)addModelViewMatrix:(std::pair<glm::mat4, glm::mat4>)matrix atMeshIndex:(int)meshIndex;
- (void)removeModelViewMatrixAtMeshIndex: (int)meshIndex atViewIndex:(int)viewIndex;
- (std::pair<glm::mat4, glm::mat4>)getModelViewMatrix:(int)meshIndex atViewIndex:(int)viewIndex;

- (BOOL)hasValidCapture:(int)meshIndex;
- (int)numViews:(int)meshIndex;

- (void)encodeWithCoder:(NSCoder *)aCoder;
- (id)copyWithZone:(NSZone *)zone;

@end
