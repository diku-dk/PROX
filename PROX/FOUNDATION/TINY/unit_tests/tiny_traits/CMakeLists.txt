add_executable(
  unit_tiny_traits
  tiny_traits.cpp
  )

target_link_libraries(
  unit_tiny_traits
  tiny
  Boost::unit_test_framework
  )

add_test(
  unit_tiny_traits
  unit_tiny_traits
  )
