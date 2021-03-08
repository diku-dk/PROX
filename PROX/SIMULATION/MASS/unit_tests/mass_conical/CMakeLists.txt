add_executable(
  unit_mass_conical
  mass_conical.cpp
  )

target_link_libraries(
  unit_mass_conical 
  mass
  Boost::unit_test_framework
  )

add_test(
  unit_mass_conical
  unit_mass_conical
  )
