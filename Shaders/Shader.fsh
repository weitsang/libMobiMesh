//
//  Shader.fsh
//  OpenGLES2Renderer
//
//  Created by Ong Yuh Shin on 4/7/11.
//  Copyright 2011 National University of Singapore. All rights reserved.
//

varying mediump vec3 destinationNormal;
varying highp vec4 destinationColor;

uniform int renderWithColor;

highp vec3 computeLighting(highp vec3 lightPos, highp float shininess,
                           highp vec3 ambientMaterial, highp vec3 specularMaterial,
                           highp vec3 diffuseMaterial)
{
    highp vec3 N = normalize(destinationNormal);
    highp vec3 L = normalize(lightPos);
    highp vec3 E = vec3(0.0, 0.0, 1.0);
    highp vec3 H = normalize(L + E); // half vector for specular lighting
    highp float diffuse = max(0.0, dot(N, L));
    highp float specular = max(0.0, dot(N, H));
    specular = pow(specular, shininess);

    /** Uncomment to enable toon shading
     if (diffuse < 0.1)
        diffuse = 0.0;
     else if (diffuse < 0.3)
        diffuse = 0.3;
     else if (diffuse < 0.6)
        diffuse = 0.6;
     else
        diffuse = 1.0;

     specular = step(0.5, specular); */
    
    return ambientMaterial + diffuse * diffuseMaterial + specular * specularMaterial;
}

void main()
{
    if (renderWithColor == 1)
    {
        gl_FragColor = destinationColor;
    }
    else
    {
        highp vec3 lightPos = vec3(0.0, 0.0, 1.0);
        highp float shininess = 25.6;
        highp vec3 ambientMaterial = vec3(0.2125, 0.1275, 0.054);
        highp vec3 diffuseMaterial = vec3(0.714, 0.4284, 0.18144);
        highp vec3 specularMaterial = vec3(0.393548, 0.271906, 0.166721);
        gl_FragColor = vec4(computeLighting(lightPos, shininess,
                                            ambientMaterial, specularMaterial, diffuseMaterial),
                            1.0);
    }
}
