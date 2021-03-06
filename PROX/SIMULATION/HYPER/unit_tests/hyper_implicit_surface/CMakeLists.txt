INCLUDE_DIRECTORIES(
  ${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/UTIL/UTIL/include 
  ${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/TINY/TINY/include 
  ${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/SPARSE/SPARSE/include 
  ${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/GEOMETRY/GEOMETRY/include 
  ${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/MESH_ARRAY/MESH_ARRAY/include 
  ${PROJECT_SOURCE_DIR}/PROX/SIMULATION/KDOP/KDOP/include 
  ${PROJECT_SOURCE_DIR}/PROX/SIMULATION/HYPER/HYPER/include 
  ${Boost_INCLUDE_DIRS}
  )

IF(HAS_DIKUCL)
  INCLUDE_DIRECTORIES(  ${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/DIKUCL/DIKUCL/include 
)
ENDIF(HAS_DIKUCL)

ADD_EXECUTABLE(
  unit_hyper_implicit_surface
  hyper_implicit_surface.cpp
  )

TARGET_LINK_LIBRARIES(
  unit_hyper_implicit_surface
  util
  tiny
  sparse
  geometry
  mesh_array
  hyper
  kdop
  ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
  )

IF(HAS_DIKUCL)
  TARGET_LINK_LIBRARIES(unit_hyper_implicit_surface dikucl)
ENDIF()

ADD_TEST(
  unit_hyper_implicit_surface
  unit_hyper_implicit_surface
  )

