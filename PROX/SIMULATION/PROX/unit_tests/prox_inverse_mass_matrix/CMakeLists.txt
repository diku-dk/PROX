INCLUDE_DIRECTORIES( 
  ${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/UTIL/UTIL/include 
  ${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/TINY/TINY/include 
  ${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/SPARSE/SPARSE/include 
  ${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/GEOMETRY/GEOMETRY/include 
  ${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/MESH_ARRAY/MESH_ARRAY/include
  ${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/CONVEX/CONVEX/include 
  ${PROJECT_SOURCE_DIR}/PROX/SIMULATION/NARROW/NARROW/include
  ${PROJECT_SOURCE_DIR}/PROX/SIMULATION/KDOP/KDOP/include
  ${PROJECT_SOURCE_DIR}/PROX/SIMULATION/BROAD/BROAD/include 
  ${PROJECT_SOURCE_DIR}/PROX/SIMULATION/PROX/PROX/include 
  ${Boost_INCLUDE_DIRS} 
)

IF(HAS_DIKUCL)
  INCLUDE_DIRECTORIES(
      ${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/DIKUCL/DIKUCL/include 
    )
ENDIF(HAS_DIKUCL)

add_executable(
  unit_prox_inverse_mass_matrix
  prox_inverse_mass_matrix.cpp
  )

target_link_libraries(
  unit_prox_inverse_mass_matrix
  util
  tiny
  sparse
  geometry
  convex
  mesh_array
  broad
  narrow
  kdop
  prox
  Boost::unit_test_framework
)

IF(HAS_DIKUCL)
  target_link_libraries(unit_prox_inverse_mass_matrix dikucl)
ENDIF()

add_test(
  unit_prox_inverse_mass_matrix
  unit_prox_inverse_mass_matrix
  )
