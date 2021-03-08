add_executable(
  unit_prox_mass_block
  prox_mass_block.cpp
  )

target_link_libraries(
  unit_prox_mass_block
  prox
  Boost::unit_test_framework
  )

add_test(
  unit_prox_mass_block
  unit_prox_mass_block
  )
