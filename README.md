# PROX
Research related to computational modeling on projects for contact, biomechanics, soft robots and more

## Install Tools

Make sure you got latest version of CMake installed on your system and your favorite IDE/compiler.

## Install Dependencies

Before compiling and linking one needs to install Boost, Glew, Qt, and Atlas.

Atlas is not a strong dependency and could be omitted. However, it is
there then PROX/FOUNDATION/BIG will take advantage of it.

### Dependencies on Apple

Atlas is already bundled with Mac OS, so no need to do anything for this.

Other dependencies can be easy installed using macport (you need to install macport if you do not already have it)

    sudo port install boost@1.53.0_1+no_single+no_static+python27+universal
    sudo port install glew@1.9.0_0+universal


### Dependencies on Linux

#### Ubuntu based versions

Install dependencies:

    $ sudo apt-get install build-essential cmake cmake-curses-gui libboost-all-dev libglew-dev libxi-dev libxmu-dev --assume-yes
    $ sudo apt-get install libatlas-dev libatlas3gf-base libatlas-base-dev --assume-yes

This has been tested on a Precise64 Virtual Box and results in the following versions:

- gcc (Ubuntu/Linaro 4.6.3-1ubuntu5)
- cmake version 2.8.7
- libglut.so.3.9.0, libGLU.so.1.3.08004, libGLEW.so.1.6.0, libGL.so.1.2
- Boost 1.46.1
- libXi.so.6.1.0, libXmu.so.6.2.0
- libatlas.so.3gf.0 with libcblas.so.3gf.0 and liblapack_atlas.so.3gf.0

#### Dependencies on Windows

Unknown as of this writing

## How to Compile and Link

CMake should be used for generating make/solution files.

* Open the cmake gui application

* Browse the source path to the top most folder containing a CMakeLists.txt file.

* For binary destination use the same folder location but add a
subfolder to contain the binary tree. You can name it whatever you
want. For instance "build". Press configure button, resolve any conflicts or issues.

* If you want to have unit-testing find the cmake variable ENABLE_UNIT_TESTS and set it to true (=1)

* Press configure button again until no red variables appear. Now press generate button.

* Goto the build folder and locate your solution/make files and use these to compile and link the sub project.


# Known Issues and Wanted Features
## ID 1

Merge SIMULATORS into PROX... Currently it only contains the
prox::Engine implementation so it seems unneeded to have a whole
separate library for this

## ID 2

Add proper ProxEngine params for gravity and damping parameters so
they too use the set_parameter interface

## ID 3

Add proper material+lighting setup to RIGID_BODY_XML2MEL application

## ID 4

Add run-time control of all input/output folders to all applications
(not sure how we want this to work?)

## ID 5

Eliminate redundant code/functionality between GUI and CMD applications

## ID 6

Eliminate redundant code/functionality between RIGID_BODY_GUI and
RIGID_BODY_GUI_CHALK applications

## ID 7

Make sure all unit-tests are working, currently these tests fails

	 11 - unit_big_gmres (Failed)
	 88 - unit_geometry_contacts_obb_cylinder (Failed)
	106 - unit_mesh_array_compute_distance_map (Failed)
	183 - unit_prox_mappings (Failed)
	184 - unit_prox_solver (Failed)

Investigations:

unit_big_gmres only fails when used with the CMake RUN_TESTS
target. When invoked individually it does not fail?

The unit_geometry_contacts_obb_cylinder is not implemented at all.

The algorithm for computing distance fields is not working, hence
unit_mesh_array_compute_distance_map fails (See bug ID 11).

The unit_prox_mappings seems to fail due too aggressive
tolerance settings for GJK?

It is not clear why unit_prox_solver fails?

## ID 8

Investigate if we can drop the CMake setup parts

CMake Warning at CMakeLists.txt:37 (MESSAGE):
  Forcing to use CMAKE_CXX_FLAGS = -Wall -Wno-unused-local-typedef because
  Kenny is lazy


CMake Warning at CMakeLists.txt:52 (MESSAGE):
  Forcing to use CMAKE_OSX_DEPLOYMENT_TARGET = 10.10 because Kenny is lazy

## ID 9

RIGID_BODY_GUI and RIGID_BODY_GUI_CHALK loads cfg-files twice. The
reason is that first load_from_config reads in the cfg file in order
to set up application logic. Later prox engine loads the cfg-files
again when trying to set its parameters from a cfg-file.

This seems a little "wasteful" to parse the same files
twice... However, it make prox engine have a clean interface with
minimal usage of data types.

## ID 10

Investigate why CMake fails to handle BLAS/LAPACK, currently we have
the error messages

A library with BLAS API not found. Please specify library location.
LAPACK requires BLAS
A library with LAPACK API not found. Please specify library location.

## ID 11

Debug mesh_array_compute_distance_map so it works
Complete unit-test of mesh_array_compute_surface_map

## ID 12

Add a function to compute area weighted face normals of a T3Mesh

## ID 13

Add functions to compute quality measures of T4Meshes

## ID 14

Add functionality to bind (couple) a T3Mesh to a T4Mesh
