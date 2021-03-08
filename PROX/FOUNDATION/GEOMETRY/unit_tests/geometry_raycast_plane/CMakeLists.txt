add_executable(
  unit_geometry_raycast_plane
  geometry_raycast_plane.cpp
  )

target_link_libraries(
  unit_geometry_raycast_plane
  geometry
  Boost::unit_test_framework
  )

add_test(
  unit_geometry_raycast_plane
  unit_geometry_raycast_plane
  )
