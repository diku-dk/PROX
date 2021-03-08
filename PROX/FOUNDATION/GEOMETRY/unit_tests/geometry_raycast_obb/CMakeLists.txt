add_executable(
  unit_geometry_raycast_obb
  geometry_raycast_obb.cpp
  )

target_link_libraries(
  unit_geometry_raycast_obb
  geometry
  Boost::unit_test_framework
  )

add_test(
  unit_geometry_raycast_obb
  unit_geometry_raycast_obb
  )
