add_executable(
  unit_geometry_inside_triangle
  geometry_inside_triangle.cpp
  )

target_link_libraries(
  unit_geometry_inside_triangle
  geometry
  Boost::unit_test_framework
  )

add_test( 
  unit_geometry_inside_triangle
  unit_geometry_inside_triangle
  )
