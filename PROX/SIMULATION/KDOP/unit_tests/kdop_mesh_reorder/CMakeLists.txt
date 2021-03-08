INCLUDE_DIRECTORIES( 
  ${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/TINY/TINY/include 
  ${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/GEOMETRY/GEOMETRY/include 
  ${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/MESH_ARRAY/MESH_ARRAY/include 
  ${PROJECT_SOURCE_DIR}/PROX/SIMULATION/KDOP/KDOP/include 
  ${Boost_INCLUDE_DIRS}
  )

add_executable(
  unit_kdop_mesh_reorder
  kdop_mesh_reorder.cpp
  )

target_link_libraries(
  unit_kdop_mesh_reorder
  tiny
  geometry
  mesh_array
  tetgen
  ${Boost_SYSTEM_LIBRARY}
  ${Boost_FILESYSTEM_LIBRARY}
  Boost::unit_test_framework
  )

add_test( 
  unit_kdop_mesh_reorder
  unit_kdop_mesh_reorder
  )

