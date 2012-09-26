//
//  Shader.vsh
//  OpenGLES2Renderer
//
//  Created by Ong Yuh Shin on 4/7/11.
//  Copyright 2011 National University of Singapore. All rights reserved.
//

attribute vec3 position;
attribute vec3 normal;
attribute vec4 color;

varying vec3 destinationNormal;
varying vec4 destinationColor;

uniform mat4 modelviewProjection;
uniform mat4 normalMatrix;

void main()
{
    destinationNormal = normalize(vec3(normalMatrix * vec4(normal, 0.0)));
    destinationColor = color;

    // Multiply matrices
    gl_Position = modelviewProjection * vec4(position, 1.0);
}
