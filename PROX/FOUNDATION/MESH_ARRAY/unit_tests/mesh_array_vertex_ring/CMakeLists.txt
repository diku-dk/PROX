add_executable(
  unit_mesh_array_vertex_ring 
  mesh_array_vertex_ring.cpp
  )

target_link_libraries(
  unit_mesh_array_vertex_ring
  mesh_array
  Boost::unit_test_framework
  )

add_test( 
  unit_mesh_array_vertex_ring
  unit_mesh_array_vertex_ring
  )
