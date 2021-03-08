add_executable(
  unit_tiny_power2
  tiny_power2.cpp
  )

target_link_libraries(
  unit_tiny_power2
  tiny
  Boost::unit_test_framework
  )

add_test(
  unit_tiny_power2
  unit_tiny_power2
  )
