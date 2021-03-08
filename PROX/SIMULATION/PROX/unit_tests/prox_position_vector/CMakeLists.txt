add_executable(
  unit_prox_position_vector
  prox_position_vector.cpp
  )

target_link_libraries(
  unit_prox_position_vector
  prox
  Boost::unit_test_framework
  )

add_test(
  unit_prox_position_vector
  unit_prox_position_vector
  )
