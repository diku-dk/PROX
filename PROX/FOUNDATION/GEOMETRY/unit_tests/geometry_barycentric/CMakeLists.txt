add_executable(
  unit_geometry_barycentric
  geometry_barycentric.cpp
  )

target_link_libraries(
  unit_geometry_barycentric
  geometry
  Boost::unit_test_framework
  )

add_test( 
  unit_geometry_barycentric
  unit_geometry_barycentric
  )
