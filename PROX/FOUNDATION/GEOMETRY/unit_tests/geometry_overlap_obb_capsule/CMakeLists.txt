add_executable(
  unit_geometry_overlap_obb_capsule
  geometry_overlap_obb_capsule.cpp
  )

target_link_libraries(
  unit_geometry_overlap_obb_capsule
  geometry
  Boost::unit_test_framework
  )

add_test( 
  unit_geometry_overlap_obb_capsule
  unit_geometry_overlap_obb_capsule
  )
