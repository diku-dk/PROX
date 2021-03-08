add_executable(
  unit_geometry_make_intersection
  geometry_make_intersection.cpp
  )

target_link_libraries(
  unit_geometry_make_intersection
  geometry
  Boost::unit_test_framework
  )

add_test( 
  unit_geometry_make_intersection
  unit_geometry_make_intersection
  )
