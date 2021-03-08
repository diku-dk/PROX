add_executable(
  unit_prox_math_policy_functions
  prox_math_policy_functions.cpp
  )

target_link_libraries(
  unit_prox_math_policy_functions
  prox
  Boost::unit_test_framework
  )

add_test(
  unit_prox_math_policy_functions
  unit_prox_math_policy_functions
  )
