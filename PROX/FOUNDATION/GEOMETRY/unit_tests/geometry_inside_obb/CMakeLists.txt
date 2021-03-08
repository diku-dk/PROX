add_executable(
  unit_geometry_inside_obb
  geometry_inside_obb.cpp
  )

target_link_libraries(
  unit_geometry_inside_obb
  geometry
  Boost::unit_test_framework
  )

add_test(
  unit_geometry_inside_obb
  unit_geometry_inside_obb
  )
