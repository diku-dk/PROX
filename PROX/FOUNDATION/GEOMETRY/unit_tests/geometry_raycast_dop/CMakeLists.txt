add_executable(
  unit_geometry_raycast_dop
  geometry_raycast_dop.cpp
  )

target_link_libraries(
  unit_geometry_raycast_dop
  geometry
  Boost::unit_test_framework
  )

add_test(
  unit_geometry_raycast_dop
  unit_geometry_raycast_dop
  )
