//
//  ModelViewCaptureContainer.mm
//  MobiMesh
//
//  Created by Ong Yuh Shin on 5/23/11.
//  Copyright 2011 National University of Singapore. All rights reserved.
//

#import "ModelViewCaptureContainer.h"


#define kVectorSize @"numMeshes"
#define kNumViews(meshIndex) [NSString stringWithFormat:@"numViews_%d", meshIndex]
#define kLookat(meshIndex, viewIndex, col, row) \
    [NSString stringWithFormat:@"lookat_%d_%d_%d_%d", meshIndex, viewIndex, col, row ]
#define kTransformation(meshIndex, viewIndex, col, row) \
    [NSString stringWithFormat:@"transformation_%d_%d_%d_%d", meshIndex, viewIndex, col, row ]

@implementation ModelViewCaptureContainer

- (id)initWithSize:(int)size
{
    self = [super init];
    modelViewMatrices.resize(size);
    return self;
}

- (id)initWithCoder:(NSCoder *)decoder
{
    if ((self = [super init]))
    {
        // Read the number of meshes
        int vectorSize = [decoder decodeIntForKey:kVectorSize];
        modelViewMatrices.resize(vectorSize);

        for (int meshIndex=0; meshIndex<vectorSize; ++meshIndex)
        {
            // Read the number of views for the current mesh
            int numViews = [decoder decodeIntForKey:kNumViews(meshIndex)];
            modelViewMatrices[meshIndex].resize(numViews);

            for (int viewIndex=0; viewIndex<numViews; ++viewIndex)
            {
                std::pair<glm::mat4, glm::mat4> modelViewMatrix;
                // Read the lookat matrix
                for (int col=0; col<4; ++col)
                    for (int row=0; row<4; ++row)
                    {
                        modelViewMatrix.first[col][row] =
                            [decoder decodeFloatForKey:kLookat(meshIndex, viewIndex, col, row)];
                    }
                // Read the transformation matrix
                for (int col=0; col<4; ++col)
                {
                    for (int row=0; row<4; ++row)
                    {
                        modelViewMatrix.second[col][row] =
                            [decoder decodeFloatForKey:kTransformation(meshIndex, viewIndex, col, row)];
                    }
                }
                // Store this view
                modelViewMatrices[meshIndex][viewIndex] = modelViewMatrix;
            }
        }
    }

    return self;
}

- (void)encodeWithCoder:(NSCoder *)aCoder
{
    // Store the number of meshes
    [aCoder encodeInt:modelViewMatrices.size() forKey:kVectorSize];
    for (int meshIndex=0; meshIndex<modelViewMatrices.size(); ++meshIndex)
    {
        // For each mesh, store the number of views
        [aCoder encodeInt:modelViewMatrices[meshIndex].size() forKey:kNumViews(meshIndex)];

        for (int viewIndex=0; viewIndex<modelViewMatrices[meshIndex].size(); ++viewIndex)
        {

            // Store the lookat matrix
            for (int col=0; col<4; ++col)
                for (int row=0; row<4; ++row)
                {
                    [aCoder encodeFloat:modelViewMatrices[meshIndex][viewIndex].first[col][row]
                                 forKey:kLookat(meshIndex, viewIndex, col, row)];
                }

            // Store the transformation matrix
            for (int col=0; col<4; ++col)
                for (int row=0; row<4; ++row)
                {
                    [aCoder encodeFloat:modelViewMatrices[meshIndex][viewIndex].second[col][row]
                                 forKey:kTransformation(meshIndex, viewIndex, col, row)];
                }
        }
    }
}

- (id)copyWithZone:(NSZone *)zone
{
    ModelViewCaptureContainer* copy = [[[self class] allocWithZone:zone] init];
    copy->modelViewMatrices = modelViewMatrices;
    return copy;
}

- (BOOL)setModelViewMatrix:(std::pair<glm::mat4, glm::mat4>)matrix atMeshIndex:(int)meshIndex atViewIndex:(int)viewIndex
{
    if (meshIndex >= modelViewMatrices.size() || index < 0 ||
        viewIndex >= modelViewMatrices[meshIndex].size() || viewIndex < 0)
        return NO; // no such matrix

    modelViewMatrices[meshIndex][viewIndex] = matrix;
    return YES;
}

- (int)addModelViewMatrix:(std::pair<glm::mat4, glm::mat4>)matrix atMeshIndex:(int)meshIndex
{
    if (meshIndex >= modelViewMatrices.size() || meshIndex < 0)
        return -1;

    modelViewMatrices[meshIndex].push_back(matrix);
    return modelViewMatrices[meshIndex].size() - 1;
}

- (void)removeModelViewMatrixAtMeshIndex: (int)meshIndex atViewIndex:(int)viewIndex
{
    if (meshIndex >= modelViewMatrices.size() || meshIndex < 0)
        return;

    std::vector<std::pair<glm::mat4, glm::mat4> >& captures = modelViewMatrices[meshIndex];
    if (viewIndex >= captures.size() || viewIndex < 0)
        return;

    captures.erase(captures.begin()+viewIndex);
}

- (std::pair<glm::mat4, glm::mat4>)getModelViewMatrix:(int)meshIndex atViewIndex:(int)viewIndex
{
    if (meshIndex >= modelViewMatrices.size() || index < 0)
        return std::pair<glm::mat4, glm::mat4>();
    else if (viewIndex >= modelViewMatrices[meshIndex].size() || viewIndex < 0)
        return std::pair<glm::mat4, glm::mat4>();

    return modelViewMatrices[meshIndex][viewIndex];
}

- (BOOL)hasValidCapture:(int)meshIndex
{
    if (meshIndex >= modelViewMatrices.size() || index < 0 || modelViewMatrices[meshIndex].size() == 0)
        return NO;
    else
        return YES;
}

- (int)numViews:(int)meshIndex
{
    if (meshIndex >= modelViewMatrices.size() || index < 0)
        return 0;
    else
        return modelViewMatrices[meshIndex].size();
}

@end
