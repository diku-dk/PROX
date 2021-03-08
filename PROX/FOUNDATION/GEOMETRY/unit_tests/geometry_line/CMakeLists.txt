add_executable(
  unit_geometry_line
  geometry_line.cpp
  )

target_link_libraries(
  unit_geometry_line
  geometry
  Boost::unit_test_framework
  )

add_test( 
  unit_geometry_line
  unit_geometry_line
  )
