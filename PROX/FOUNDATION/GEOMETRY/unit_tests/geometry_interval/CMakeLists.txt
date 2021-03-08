add_executable(
  unit_geometry_interval 
  geometry_interval.cpp
  )

target_link_libraries(
  unit_geometry_interval
  geometry
  Boost::unit_test_framework
  )

add_test( 
  unit_geometry_interval
  unit_geometry_interval
  )
