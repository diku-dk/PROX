INCLUDE_DIRECTORIES( 
  ${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/MESH_ARRAY/MESH_ARRAY/include 
  ${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/TINY/TINY/include
  ${Boost_INCLUDE_DIRS}
  )

add_executable(
  unit_mesh_array_compute_distance_map
  mesh_array_compute_distance_map.cpp
  )

target_link_libraries(
  unit_mesh_array_compute_distance_map
  mesh_array
  tiny
  tetgen
  ${Boost_SYSTEM_LIBRARY}
  ${Boost_FILESYSTEM_LIBRARY}
  Boost::unit_test_framework
  )

add_test( 
  unit_mesh_array_compute_distance_map
  unit_mesh_array_compute_distance_map
  )

