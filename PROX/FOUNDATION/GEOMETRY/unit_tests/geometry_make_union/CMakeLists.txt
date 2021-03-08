add_executable(
  unit_geometry_make_union
  geometry_make_union.cpp
  )

target_link_libraries(
  unit_geometry_make_union
  geometry
  Boost::unit_test_framework
  )

add_test( 
  unit_geometry_make_union
  unit_geometry_make_union
  )
