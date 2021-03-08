add_executable(
  unit_geometry_raycast_aabb
  geometry_raycast_aabb.cpp
  )

target_link_libraries(
  geometry
  Boost::unit_test_framework
  )

add_test(
  unit_geometry_raycast_aabb
  unit_geometry_raycast_aabb
  )
