add_executable(
  unit_mesh_array_triangle 
  mesh_array_triangle.cpp
  )

target_link_libraries(
  unit_mesh_array_triangle
  mesh_array
  Boost::unit_test_framework
  )

add_test( 
  unit_mesh_array_triangle
  unit_mesh_array_triangle
  )
