add_executable(
  unit_mass_rotate
  mass_rotate.cpp
  )

target_link_libraries(
  unit_mass_rotate
  mass
  Boost::unit_test_framework
  )

add_test(
  unit_mass_rotate
  unit_mass_rotate
  )
