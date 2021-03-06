cmake_minimum_required(VERSION 3.12)

project(PROX VERSION 0.1
                   DESCRIPTION "A header only c++ library for simulating contact between rigid or soft bodies"
                   HOMEPAGE_URL "https://github.com/diku-dk/PROX")

set(CMAKE_CONFIGURATION_TYPES "Debug;Release" CACHE STRING "" FORCE)

add_definitions(-D_USE_MATH_DEFINES)

#-----------------------------------------------------------------------------
#
# Conan support for PROX
#
# 2021-03-06 Kenny: Preparing for conan
# # Add connan dependencies, conan targets will be used.
# if(EXISTS ${CMAKE_BINARY_DIR}/conanbuildinfo.cmake)
#  include(${CMAKE_BINARY_DIR}/conanbuildinfo.cmake)
#  conan_basic_setup(TARGETS NO_OUTPUT_DIRS)
#  set(Boost_USE_STATIC_LIBS ON)
# endif()

#-----------------------------------------------------------------------------
#
# Avoid in-source builds.
#
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_BINARY_DIR)
  message(
    FATAL_ERROR
    "CMake generation for PROX is not allowed within the source directory!"
    )
endif()

#-----------------------------------------------------------------------------
#
# Global variables that control the behaviour of CMake. Use these in
# the GUI to turn on/off different kind of things.
#
option(PROX_ENABLE_PROFILING "Build demos" OFF)
option(PROX_ENABLE_UNIT_TESTS "Build unit test" OFF)

if(PROX_ENABLE_PROFILING)
  message("Profiling is...........................ON")
  add_definitions(-DUSE_PROFILING)
else()
  message("Profiling is...........................OFF")
endif()

if(PROX_ENABLE_UNIT_TESTS)
  message("Unit-tests are.........................ON")
  enable_testing()
else()
  message("Unit-tests are.........................OFF")
endif()

#-----------------------------------------------------------------------------
#
# Find all libraries used by PROX
#
SET(Boost_USE_STATIC_LIBS ON)

FIND_PACKAGE(
  Boost
  COMPONENTS
  unit_test_framework
  system
  filesystem
  REQUIRED
  )

FIND_PACKAGE(OpenGL REQUIRED)
FIND_PACKAGE(GLEW REQUIRED)

SET( TMP_CMAKE_PREFIX_PATH  ${CMAKE_PREFIX_PATH} )
SET( CMAKE_PREFIX_PATH "/Users/kenny/Qt/5.3/clang_64/" ${TMP_CMAKE_PREFIX_PATH} )
SET( CMAKE_PREFIX_PATH "/Users/kennyerleben/Qt/5.3/clang_64/" ${TMP_CMAKE_PREFIX_PATH} )
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




ADD_SUBDIRECTORY(PROX)
