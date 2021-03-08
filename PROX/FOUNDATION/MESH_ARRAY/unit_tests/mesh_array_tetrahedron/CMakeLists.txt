add_executable(
  unit_mesh_array_tetrahedron 
  mesh_array_tetrahedron.cpp
  )

target_link_libraries(
  unit_mesh_array_tetrahedron
  mesh_array
  Boost::unit_test_framework
  )

add_test( 
  unit_mesh_array_tetrahedron
  unit_mesh_array_tetrahedron
  )
