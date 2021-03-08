add_executable(
  unit_mass_density
  mass_density.cpp
  )

target_link_libraries(
  unit_mass_density
  mass
  Boost::unit_test_framework
  )

add_test(
  unit_mass_density
  unit_mass_density
  )
