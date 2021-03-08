add_executable(
  unit_geometry_ray
  geometry_ray.cpp
  )

target_link_libraries(
  unit_geometry_ray
  geometry
  Boost::unit_test_framework
  )

add_test( 
  unit_geometry_ray
  unit_geometry_ray
  )
