add_executable(
  unit_geometry_contacts_obb_sphere
  geometry_contacts_obb_sphere.cpp
  )

target_link_libraries(
  unit_geometry_contacts_obb_sphere
  geometry
  Boost::unit_test_framework
  )

add_test( 
  unit_geometry_contacts_obb_sphere
  unit_geometry_contacts_obb_sphere
  )
