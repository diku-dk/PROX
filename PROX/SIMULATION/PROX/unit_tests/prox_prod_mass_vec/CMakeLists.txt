add_executable(
  unit_prox_prod_mass_vec
  prox_prod_mass_vec.cpp
  )

target_link_libraries(
  unit_prox_prod_mass_vec 
  prox
  Boost::unit_test_framework
  )

add_test(
  unit_prox_prod_mass_vec
  unit_prox_prod_mass_vec
  )
