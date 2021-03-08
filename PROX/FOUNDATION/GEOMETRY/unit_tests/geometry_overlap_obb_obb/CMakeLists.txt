add_executable(
  unit_geometry_overlap_obb_obb
  geometry_overlap_obb_obb.cpp
  )

target_link_libraries(
  unit_geometry_overlap_obb_obb
  geometry
  Boost::unit_test_framework
  )

add_test( 
  unit_geometry_overlap_obb_obb
  unit_geometry_overlap_obb_obb
  )
