add_executable(
  unit_prox_prod_mass_jac
  prox_prod_mass_jac.cpp
  )

target_link_libraries(
  unit_prox_prod_mass_jac
  prox
  Boost::unit_test_framework
  )

add_test(
  unit_prox_prod_mass_jac
  unit_prox_prod_mass_jac
  )
