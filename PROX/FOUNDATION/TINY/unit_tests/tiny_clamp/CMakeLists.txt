add_executable(
  unit_tiny_clamp
  tiny_clamp.cpp
  )

target_link_libraries(
  unit_tiny_clamp
  tiny
  Boost::unit_test_framework
  )

add_test(
  unit_tiny_clamp
  unit_tiny_clamp
  )
