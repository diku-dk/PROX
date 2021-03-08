add_executable(
  unit_mesh_array_adjacency 
  mesh_array_adjacency.cpp
  )

target_link_libraries(
  unit_mesh_array_adjacency
  mesh_array
  Boost::unit_test_framework
  )

add_test( 
  unit_mesh_array_adjacency
  unit_mesh_array_adjacency
  )
