This file contains a description of the directory structure, as well as
information on how to build the library and application(s).

The library and projects are set up on xCode with the project
MobiMesh.xcodeproj.

Directory Structure
-------------------
The following are the list of directories and what they contain:

Classes
  Files for compiling and running MobiMesh application. These are class
  and interface files common to running on both the iPad and iPod Touch.

MeshFormatConverters
  Contains the code for converting from traditional .pm files to .vpm,
  .tvpm and .tvpm64 formats.

MobiMeshIPad
  Files for compiling and running MobiMesh application specifically for
  iPad. This includes the main.m and the nib files.

Shaders
  OpenGL ES 2.0 shaders used by the MobiMesh application as well as the
  sample tutorial code.

TechReport
  Latex source files and relevant diagrams in .ps format for compiling
  the technical report.

Tutorial
  Simple tutorial on how to use libMobiMesh.

TVPMRenderer
  Files for the tutorial's sample code, including interface/class files
  and nib files.

glm
  Contains the GLM library headers required for building the MobiMesh
  library.

libMobiMesh
  Source code fo the MobiMesh library.


Building libMobiMesh
--------------------
1) Open MobiMesh.xcodeproj
2) Select the libMobiMesh scheme/target
3) Click run

Building and Running MobiMesh Application
-----------------------------------------
1) Open MobiMesh.xcodeproj
2) If building for iPad, select MobiMeshIPad scheme/target
   If building for iPod Touch/iPhone, select MobiMesh scheme/target
3) Click run

Building and Running Sample TVPMRenderer
----------------------------------------
1) Open MobiMesh.xcodeproj
2) Select the TVPMRenderer scheme/target
3) Click run

MobiMesh.xcodeproj Organization
-------------------------------
The xCode project MobiMesh contains 4 targets:
  - MobiMesh: iPod Touch/iPhone application
  - libMobiMesh: the MobiMesh C++ library
  - MobiMeshIPad: iPad application
  - TVPMRenderer: Sample code for rendering .tvpm meshes

The project is organized into the following groups:
  - glm: Source code/headers for GLM library
  - libMobiMesh: Source code/headers for libMobiMesh
    - Exceptions: Exception classes for libMobiMesh
  - Classes: Source code/headers for MobiMesh target
  - Shaders: Vertex and fragment shader for OpenGL ES 2.0
  - Other Sources: Misc files (xCode generated) for MobiMesh target
  - Resources: Mesh files, icon and nib files for MobiMesh target
  - MobiMeshIPad: Nib files for MobiMeshIPad target plus other
    supporting files
  - TVPMRenderer: Source code/headers, nib and supporting files for
    TVPMRenderer target.
