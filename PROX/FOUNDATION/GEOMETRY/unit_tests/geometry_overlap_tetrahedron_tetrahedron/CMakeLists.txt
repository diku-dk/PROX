add_executable(
  unit_geometry_overlap_tetrahedron_tetrahedron
  geometry_overlap_tetrahedron_tetrahedron.cpp
  )

target_link_libraries(
  unit_geometry_overlap_tetrahedron_tetrahedron
  geometry
  Boost::unit_test_framework
  )

add_test( 
  unit_geometry_overlap_tetrahedron_tetrahedron
  unit_geometry_overlap_tetrahedron_tetrahedron
  )
