add_executable(
  unit_mesh_array_mesh 
  mesh_array_mesh.cpp
  )

target_link_libraries(
  unit_mesh_array_mesh
  mesh_array
  Boost::unit_test_framework
  )

add_test( 
  unit_mesh_array_mesh
  unit_mesh_array_mesh
  )
