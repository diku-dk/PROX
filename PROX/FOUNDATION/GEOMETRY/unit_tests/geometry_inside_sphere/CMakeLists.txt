add_executable(
  unit_geometry_inside_sphere
  geometry_inside_sphere.cpp
  )

target_link_libraries(
  unit_geometry_inside_sphere
  geometry
  Boost::unit_test_framework
  )

add_test(
  unit_geometry_inside_sphere
  unit_geometry_inside_sphere
  )
