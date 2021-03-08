add_executable(
  unit_tiny_basic_math_types
  tiny_basic_math_types.cpp
  )

target_link_libraries(
  unit_tiny_basic_math_types
  tiny
  Boost::unit_test_framework
  )

add_test(
  unit_tiny_basic_math_types
  unit_tiny_basic_math_types
  )
