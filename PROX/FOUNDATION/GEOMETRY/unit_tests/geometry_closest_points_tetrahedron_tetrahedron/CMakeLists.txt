add_executable(
  unit_geometry_closest_points_tetrahedron_tetrahedron
  geometry_closest_points_tetrahedron_tetrahedron.cpp
  )

target_link_libraries(
  unit_geometry_closest_points_tetrahedron_tetrahedron
  geometry
  Boost::unit_test_framework
  )

add_test( 
  unit_geometry_closest_points_tetrahedron_tetrahedron
  unit_geometry_closest_points_tetrahedron_tetrahedron
  )
