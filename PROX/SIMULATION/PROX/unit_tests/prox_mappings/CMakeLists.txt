add_executable(
  unit_prox_mappings
  prox_mappings.cpp
  )

target_link_libraries(
  unit_prox_mappings
  prox
  Boost::unit_test_framework
  )

add_test(
  unit_prox_mappings
  unit_prox_mappings
  )
