add_executable(
  unit_mass_sphere
  mass_sphere.cpp
  )

target_link_libraries(
  unit_mass_sphere 
  mass
  Boost::unit_test_framework
  )

add_test(
  unit_mass_sphere
  unit_mass_sphere
  )
