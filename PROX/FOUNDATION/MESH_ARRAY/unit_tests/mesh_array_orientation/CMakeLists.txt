add_executable(
  unit_mesh_array_orientation
  mesh_array_orientation.cpp
  )

target_link_libraries(
  unit_mesh_array_orientation
  Boost::unit_test_framework
  mesh_array
  )

add_test( 
  unit_mesh_array_orientation
  unit_mesh_array_orientation
  )
