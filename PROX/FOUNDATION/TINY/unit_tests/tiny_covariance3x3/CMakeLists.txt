add_executable(
  unit_tiny_covariance3x3
  tiny_covariance3x3.cpp
  )

target_link_libraries(
  unit_tiny_covariance3x3
  tiny
  Boost::unit_test_framework
  )

add_test(
  unit_tiny_covariance3x3
  unit_tiny_covariance3x3
  )
