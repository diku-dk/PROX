IF(HAS_DIKUCL)
  CONFIGURE_FILE(
    "${PROJECT_SOURCE_DIR}/PROX/SIMULATION/HYPER/HYPER/include/cl/hyper_cl_kernels_path.h.in"
    "${PROJECT_SOURCE_DIR}/PROX/SIMULATION/HYPER/HYPER/include/cl/hyper_cl_kernels_path.h"
    @ONLY
    )

  INCLUDE_DIRECTORIES(
      ${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/DIKUCL/DIKUCL/include 
    )
  
ENDIF(HAS_DIKUCL)

INCLUDE_DIRECTORIES( 
  ${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/TINY/TINY/include 
  ${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/SPARSE/SPARSE/include 
  ${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/MESH_ARRAY/MESH_ARRAY/include 
  ${PROJECT_SOURCE_DIR}/PROX/SIMULATION/HYPER/HYPER/include 
  ${Boost_INCLUDE_DIRS}
)

FILE(
  GLOB
  HYPER_SOURCES
  ${PROJECT_SOURCE_DIR}/PROX/SIMULATION/HYPER/HYPER/src/*.cpp
  )

FILE(
  GLOB
  HYPER_HEADERS ${PROJECT_SOURCE_DIR}/PROX/SIMULATION/HYPER/HYPER/include/*.h
  )

ADD_LIBRARY(
  hyper
  ${HYPER_SOURCES}
  ${HYPER_HEADERS}
  )

TARGET_LINK_LIBRARIES(
  hyper
  tiny
  sparse
  mesh_array
  )

IF(HAS_DIKUCL)
  TARGET_LINK_LIBRARIES(hyper dikucl)
ENDIF(HAS_DIKUCL)

SET_TARGET_PROPERTIES(hyper PROPERTIES DEBUG_POSTFIX "D" )

SET_TARGET_PROPERTIES(hyper PROPERTIES ARCHIVE_OUTPUT_DIRECTORY                "${PROJECT_SOURCE_DIR}/lib/")
SET_TARGET_PROPERTIES(hyper PROPERTIES ARCHIVE_OUTPUT_DIRECTORY_DEBUG          "${PROJECT_SOURCE_DIR}/lib/")
SET_TARGET_PROPERTIES(hyper PROPERTIES ARCHIVE_OUTPUT_DIRECTORY_RELEASE        "${PROJECT_SOURCE_DIR}/lib/")
SET_TARGET_PROPERTIES(hyper PROPERTIES ARCHIVE_OUTPUT_DIRECTORY_MINSIZEREL     "${PROJECT_SOURCE_DIR}/lib/")
SET_TARGET_PROPERTIES(hyper PROPERTIES ARCHIVE_OUTPUT_DIRECTORY_RELWITHDEBINFO "${PROJECT_SOURCE_DIR}/lib/")
