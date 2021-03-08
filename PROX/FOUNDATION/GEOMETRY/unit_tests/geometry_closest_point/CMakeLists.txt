add_executable(
  unit_geometry_closest_point
  geometry_closest_point.cpp
  )

target_link_libraries(
  unit_geometry_closest_point
  geometry
  Boost::unit_test_framework
  )

add_test( 
  unit_geometry_closest_point
  unit_geometry_closest_point
  )
