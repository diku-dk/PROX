add_executable(
  unit_mass_translate
  mass_translate.cpp
  )

target_link_libraries(
  unit_mass_translate
  mass
  Boost::unit_test_framework
  )

add_test(
  unit_mass_translate
  unit_mass_translate
  )
