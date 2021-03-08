add_executable(
  unit_geometry_raycast_tetrahedron
  geometry_raycast_tetrahedron.cpp
  )

target_link_libraries(
  unit_geometry_raycast_tetrahedron
  geometry
  Boost::unit_test_framework
  )

add_test(
  unit_geometry_raycast_tetrahedron
  unit_geometry_raycast_tetrahedron
  )

