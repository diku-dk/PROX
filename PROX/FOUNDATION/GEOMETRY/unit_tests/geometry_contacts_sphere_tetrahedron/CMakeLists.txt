add_executable(
  unit_geometry_contacts_sphere_tetrahedron
  geometry_contacts_sphere_tetrahedron.cpp
  )

target_link_libraries(
  unit_geometry_contacts_sphere_tetrahedron
  geometry
  Boost::unit_test_framework
  )

add_test( 
  unit_geometry_contacts_sphere_tetrahedron
  unit_geometry_contacts_sphere_tetrahedron
  )
