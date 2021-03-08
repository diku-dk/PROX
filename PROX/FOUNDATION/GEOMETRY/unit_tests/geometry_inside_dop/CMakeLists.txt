add_executable(
  unit_geometry_inside_dop
  geometry_inside_dop.cpp
  )

target_link_libraries(
  unit_geometry_inside_dop
  geometry
  Boost::unit_test_framework
  )

add_test(
  unit_geometry_inside_dop
  unit_geometry_inside_dop
  )
