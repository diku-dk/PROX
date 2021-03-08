add_executable(
  unit_geometry_contacts_obb_cylinder
  geometry_contacts_obb_cylinder.cpp
  )

target_link_libraries(
  unit_geometry_contacts_obb_cylinder
  geometry
  Boost::unit_test_framework
  )

add_test( 
  unit_geometry_contacts_obb_cylinder
  unit_geometry_contacts_obb_cylinder
  )
