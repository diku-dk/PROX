add_executable(
  unit_mass_capsule
  mass_capsule.cpp
  )

target_link_libraries(
  unit_mass_capsule 
  mass
  Boost::unit_test_framework
  )

add_test(
  unit_mass_capsule
  unit_mass_capsule
  )
