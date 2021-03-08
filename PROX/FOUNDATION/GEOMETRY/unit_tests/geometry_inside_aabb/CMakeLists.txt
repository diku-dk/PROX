add_executable(
  unit_geometry_inside_aabb
  geometry_inside_aabb.cpp
  )

target_link_libraries(
  unit_geometry_inside_aabb
  geometry
  Boost::unit_test_framework
  )

add_test( 
  unit_geometry_inside_aabb
  unit_geometry_inside_aabb
  )
