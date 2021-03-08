add_executable(
  unit_mass_cylinder
  mass_cylinder.cpp
  )

target_link_libraries(
  unit_mass_cylinder
  mass
  Boost::unit_test_framework
  )

add_test(
  unit_mass_cylinder
  unit_mass_cylinder
  )
