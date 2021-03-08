add_executable(
  unit_tiny_prime_numbers
  tiny_prime_numbers.cpp
  )

target_link_libraries(
  unit_tiny_prime_numbers
  tiny
  Boost::unit_test_framework
  )

add_test(
  unit_tiny_prime_numbers
  unit_tiny_prime_numbers
  )
