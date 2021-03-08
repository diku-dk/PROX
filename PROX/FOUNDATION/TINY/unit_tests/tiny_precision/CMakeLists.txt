add_executable(
  unit_tiny_precision
  tiny_precision.cpp
  )

target_link_libraries(
  unit_tiny_precision
  tiny
  Boost::unit_test_framework
  )

add_test(
  unit_tiny_precision
  unit_tiny_precision
  )
