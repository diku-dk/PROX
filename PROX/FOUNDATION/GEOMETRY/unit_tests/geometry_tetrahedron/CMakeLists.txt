add_executable(
  unit_geometry_tetrahedron
  geometry_tetrahedron.cpp
  )

target_link_libraries(
  unit_geometry_tetrahedron
  geometry
  Boost::unit_test_framework
  )

add_test( 
  unit_geometry_tetrahedron
  unit_geometry_tetrahedron
  )
