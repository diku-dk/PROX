add_executable(
  unit_geometry_plane
  geometry_plane.cpp
  )

target_link_libraries(
  unit_geometry_plane
  geometry
  Boost::unit_test_framework
  )

add_test( 
  unit_geometry_plane
  unit_geometry_plane
  )
