add_executable(
  unit_geometry_convert 
  geometry_convert.cpp
  )

target_link_libraries(
  unit_geometry_convert
  geometry
  Boost::unit_test_framework
  )

add_test( 
  unit_geometry_convert
  unit_geometry_convert
  )
