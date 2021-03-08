add_executable(
  unit_mass_box
  mass_box.cpp
  )

target_link_libraries(
  unit_mass_box 
  mass
  Boost::unit_test_framework
  )

add_test(
  unit_mass_box
  unit_mass_box
  )
