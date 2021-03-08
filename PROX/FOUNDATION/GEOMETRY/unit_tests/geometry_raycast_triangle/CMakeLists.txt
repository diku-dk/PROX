add_executable(
  unit_geometry_raycast_triangle
  geometry_raycast_triangle.cpp
  )

target_link_libraries(
  unit_geometry_raycast_triangle
  geometry
  Boost::unit_test_framework
  )

add_test(
  unit_geometry_raycast_triangle
  unit_geometry_raycast_triangle
  )
