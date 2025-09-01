#add_compile_options(-Wall -Wextra -Wpedantic -Werror -Wno-error=unused-parameter -Wno-error=unused-result -Wno-error=unused-local-typedefs -Wno-error=unused-but-set-variable -Wno-error=deprecated-copy -Wno-error=unused-function -Wno-error=unused-variable -Wno-error=deprecated-declarations -fno-omit-frame-pointer -fsanitize=address)
add_compile_options(-fno-omit-frame-pointer -fsanitize=address)
add_link_options(-fno-omit-frame-pointer -fsanitize=address)

IF(COMMAND cmake_policy)
  CMAKE_POLICY(SET CMP0003 NEW)
  CMAKE_POLICY(SET CMP0043 NEW) # QT stuff needs this?
ENDIF(COMMAND cmake_policy)

PROJECT(PROX)

find_package(Eigen3 REQUIRED)
include_directories(EIGEN3_INCLUDE_DIR)
include_directories(${CMAKE_SOURCE_DIR}/PROX/3RDPARTY/libigl/include)
##########################################################
##########################################################
##########################################################

CMAKE_MINIMUM_REQUIRED(VERSION 3.0.0 FATAL_ERROR)

##########################################################
##########################################################
##########################################################

SET(CMAKE_CONFIGURATION_TYPES "Debug;Release" CACHE STRING "" FORCE)
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -ffloat-store")

##########################################################
##########################################################
##########################################################

IF (MSVC)

  ADD_DEFINITIONS(-D_SCL_SECURE_NO_WARNINGS)
  ADD_DEFINITIONS(-D_CRT_SECURE_NO_WARNINGS)

ENDIF(MSVC)

IF(CMAKE_GENERATOR MATCHES Xcode)

  # Add platform specific stuff here!!!

ENDIF(CMAKE_GENERATOR MATCHES Xcode)

##########################################################
##########################################################
##########################################################

#--- Make sure compiler knows we want to use c++11 features if they exist     ----
#--- This works fine with "system" compiler, but NVCC will not work with this  ---

set (CMAKE_CXX_STANDARD 20)


##########################################################
##########################################################
##########################################################

SET(
  CMAKE_MODULE_PATH
  ${PROJECT_SOURCE_DIR}/CMake
  ${PROJECT_SOURCE_DIR}/cmake
  ${CMAKE_MODULE_PATH}
  )

##########################################################
##########################################################
##########################################################

SET(Boost_USE_STATIC_LIBS ON)
IF (MSVC)
  SET( BOOST_LIBRARYDIR "${Boost_INCLUDE_DIR}/lib")
ENDIF(MSVC)

FIND_PACKAGE(
  Boost
  COMPONENTS
  unit_test_framework
  system
  filesystem
  REQUIRED
  )

##########################################################
##########################################################
##########################################################

#ADD_DEFINITIONS(-DUSE_ATLAS)  #2012-09-28 Kenny: Comment/Outcomment depending on whether one wants ATLAS "acceleration"

##########################################################
##########################################################
##########################################################

ADD_DEFINITIONS(-D_USE_MATH_DEFINES)

##########################################################
##########################################################
##########################################################

IF(APPLE)

  SET(
    BLA_VENDOR
    "ATLAS"
    CACHE STRING
    "Specify what BLA vendor to use"
    )

ELSEIF(WIN32)

  # 2012-09-28 Kenny: Could be changed to inter10_32 or intel10_64lp, see documentation for FindLAPACK.cmake
  SET(
    BLA_VENDOR
    "intel"
    CACHE STRING
    "Specify what BLA vendor to use"
    )

ELSE()

   # 2013-11-23 Robert: Linux with atlas libraries installed, see README file
   SET(
     BLA_VENCOR
     "ATLAS"
     CACHE STRING
     "Specify what BLA vendor to use"
     )

ENDIF()

FIND_PACKAGE(LAPACK)

##########################################################
##########################################################
##########################################################

IF(APPLE)

  SET(
    BIG_BLAS_LIBRARY
    ${BLAS_atlas_LIBRARY}
    CACHE STRING
    "BLAS Library to use with BIG"
    )

ELSE()

  SET(
    BIG_BLAS_LIBRARY
    ${BLAS_intel_LIBRARY}
    CACHE STRING
    "BLAS Library to use with BIG"
    )

ENDIF()

##########################################################
##########################################################
##########################################################

FIND_PATH(
  Boost_BINDINGS_INCLUDE_DIR
  boost/numeric/bindings/atlas/clapack.hpp
  "${PROJECT_SOURCE_DIR}/PROX/3RDPARTY/BOOST_BINDINGS/boost_bindings/"
  DOC "What is the path where the boost bindings can be found"
  )

IF(Boost_BINDINGS_INCLUDE_DIR)

  SET( Boost_BINDINGS_FOUND 1 CACHE STRING "Set to 1 if boost bindings is found, 0 otherwise")

ELSE()

  SET( Boost_BINDINGS_FOUND 0 CACHE STRING "Set to 1 if boost bindings is found, 0 otherwise")

ENDIF()

MARK_AS_ADVANCED( Boost_BINDINGS_FOUND           )
MARK_AS_ADVANCED( Boost_BINDINGS_INCLUDE_DIR    )

##########################################################
##########################################################
##########################################################

FIND_PACKAGE(OPENCL)

IF(NOT OPENCL_FOUND)

  MESSAGE("OpenCL could not be found on the system")

ENDIF()

##########################################################
##########################################################
##########################################################

SET(
  ENABLE_DIKUCL
  1
  CACHE STRING
  "Set to 1 if DIKUCL should be enabled and 0 otherwise"
  )

IF(ENABLE_DIKUCL)

  SET(HAS_DIKUCL ${OPENCL_FOUND})

  IF(OPENCL_FOUND)

    ADD_DEFINITIONS(-DHAS_DIKUCL)

    MESSAGE("HAS_DIKUCL is turned...................ON")

  ELSE()

    MESSAGE("HAS_DIKUCL is turned...................OFF")

  ENDIF()

ENDIF()

##########################################################
##########################################################
##########################################################

SET(ENABLE_DEBUG_DIKUCL 1 CACHE STRING "Set to 1 if OpenCL debugging should be enabled and 0 otherwise")

IF(ENABLE_DIKUCL)

  IF(HAS_DIKUCL)

    IF(NOT ENABLE_DEBUG_DIKUCL)

      MESSAGE("OpenCL Debugging is turned.............OFF")
      ADD_DEFINITIONS(-DNDEBUG_DIKUCL)

    ELSE()

      MESSAGE("OpenCL Debugging is turned.............ON")

    ENDIF()

  ENDIF()

ENDIF()

##########################################################
##########################################################
##########################################################

FIND_PACKAGE(OpenGL REQUIRED)

##########################################################
##########################################################
##########################################################

FIND_PACKAGE(GLEW REQUIRED)

##########################################################
##########################################################
##########################################################

SET( TMP_CMAKE_PREFIX_PATH  ${CMAKE_PREFIX_PATH} )

SET( CMAKE_PREFIX_PATH "/Users/kenny/Qt/5.3/clang_64/" ${TMP_CMAKE_PREFIX_PATH} )
SET( CMAKE_PREFIX_PATH "/Users/kennyerleben/Qt/5.3/clang_64/" ${TMP_CMAKE_PREFIX_PATH} )
# SET( CMAKE_PREFIX_PATH "/Users/kenny/QtX/5.3/clang_64/" ${TMP_CMAKE_PREFIX_PATH} )
# SET( CMAKE_PREFIX_PATH "/home/robert/Qt/5.3/gcc_64/" ${TMP_CMAKE_PREFIX_PATH} )

FIND_PACKAGE(Qt5Widgets  REQUIRED)
FIND_PACKAGE(Qt5Core     REQUIRED)
FIND_PACKAGE(Qt5Gui      REQUIRED)
FIND_PACKAGE(Qt5OpenGL   REQUIRED)

SET( CMAKE_PREFIX_PATH ${TMP_CMAKE_PREFIX_PATH} )


##########################################################
##########################################################
##########################################################
#
# Tell that TinyXML in 3RDPARTY should use
# std::strings over char *
#
# TinyXML can work both with stl strings and char*. The
# functionality is controlled by defining the pre-processor
# directive TIXML_USE_STL or not.
#
# Note that every target that links against TinyXML
# also need to have the TIXML_USE_STL flag defined!
#
ADD_DEFINITIONS(-DTIXML_USE_STL)

##########################################################
##########################################################
##########################################################

SET(ENABLE_PROFILING     1 CACHE STRING "Set to 1 if profiling should be enabled and 0 otherwise"                )

IF(ENABLE_PROFILING)

  MESSAGE("Profiling is...........................ON")
  ADD_DEFINITIONS(-DUSE_PROFILING)

ELSE()

  MESSAGE("Profiling is...........................OFF")

ENDIF()

##########################################################
##########################################################
##########################################################

SET(ENABLE_UNIT_TESTS 1 CACHE STRING "Set to 1 if unit tests should be added to project files and 0 otherwise")

IF(ENABLE_UNIT_TESTS)

  MESSAGE("Unit-tests are.........................ON")

  ENABLE_TESTING()

ELSE()

  MESSAGE("Unit-tests are.........................OFF")

ENDIF()

##########################################################
##########################################################
##########################################################

ADD_SUBDIRECTORY(PROX)

##########################################################
##########################################################
##########################################################


