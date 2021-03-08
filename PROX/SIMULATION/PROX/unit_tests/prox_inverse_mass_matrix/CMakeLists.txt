add_executable(
  unit_prox_inverse_mass_matrix
  prox_inverse_mass_matrix.cpp
  )

target_link_libraries(
  unit_prox_inverse_mass_matrix
  prox
  Boost::unit_test_framework
  )

add_test(
  unit_prox_inverse_mass_matrix
  unit_prox_inverse_mass_matrix
  )
