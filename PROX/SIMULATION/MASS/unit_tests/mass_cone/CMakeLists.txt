add_executable(
  unit_mass_cone
  mass_cone.cpp
  )

target_link_libraries(
  unit_mass_cone
  mass
  Boost::unit_test_framework
  )

add_test(
  unit_mass_cone
  unit_mass_cone
  )
