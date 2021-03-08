add_executable(
  unit_geometry_transform
  geometry_transform.cpp
  )

target_link_libraries(
  unit_geometry_transform
  geometry
  Boost::unit_test_framework
  )

add_test( 
  unit_geometry_transform
  unit_geometry_transform
  )

