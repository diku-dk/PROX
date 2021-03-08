add_executable(
  unit_tiny_polar_decomposition
  tiny_polar_decomposition.cpp
  )

target_link_libraries(
  unit_tiny_polar_decomposition
  tiny
  Boost::unit_test_framework
  )

add_test(
  unit_tiny_polar_decomposition
  unit_tiny_polar_decomposition
  )
