add_executable(
  unit_geometry_aabb
  geometry_aabb.cpp
  )

target_link_libraries(
  unit_geometry_aabb
  geometry
  Boost::unit_test_framework
  )

add_test( 
  unit_geometry_aabb
  unit_geometry_aabb
  )
