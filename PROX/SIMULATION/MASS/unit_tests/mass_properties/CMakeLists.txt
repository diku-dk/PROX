add_executable(
  unit_mass_properties
  mass_properties.cpp
  )

target_link_libraries(
  unit_mass_properties
  mass
  Boost::unit_test_framework
  )

add_test(
  unit_mass_properties
  unit_mass_properties
  )
