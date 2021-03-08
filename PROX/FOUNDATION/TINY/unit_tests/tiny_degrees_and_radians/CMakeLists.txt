add_executable(
  unit_tiny_degrees_and_radians
  tiny_degrees_and_radians.cpp
  )

target_link_libraries(
  unit_tiny_degrees_and_radians
  tiny
  Boost::unit_test_framework
  )

add_test(
  unit_tiny_degrees_and_radians
  unit_tiny_degrees_and_radians
  )
