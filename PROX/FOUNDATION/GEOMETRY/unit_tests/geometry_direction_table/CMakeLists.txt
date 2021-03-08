add_executable(
  unit_geometry_direction_table 
  geometry_direction_table.cpp
  )

target_link_libraries(
  unit_geometry_direction_table
  geometry
  Boost::unit_test_framework
  )

add_test( 
  unit_geometry_direction_table 
  unit_geometry_direction_table 
  )
