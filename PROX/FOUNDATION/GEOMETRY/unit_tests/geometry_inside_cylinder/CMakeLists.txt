add_executable(
  unit_geometry_inside_cylinder
  geometry_inside_cylinder.cpp
  )

target_link_libraries(
  unit_geometry_inside_cylinder
  geometry
  Boost::unit_test_framework
  )

add_test(
  unit_geometry_inside_cylinder
  unit_geometry_inside_cylinder
  )
