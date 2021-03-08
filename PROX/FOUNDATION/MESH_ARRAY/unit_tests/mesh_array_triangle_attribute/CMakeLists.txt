add_executable(
  unit_mesh_array_triangle_attribute
  mesh_array_triangle_attribute.cpp
  )

target_link_libraries(
  unit_mesh_array_triangle_attribute
  mesh_array
  Boost::unit_test_framework
  )

add_test(
  unit_mesh_array_triangle_attribute
  unit_mesh_array_triangle_attribute
  )
