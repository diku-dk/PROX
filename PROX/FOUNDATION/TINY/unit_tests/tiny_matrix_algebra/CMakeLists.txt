add_executable(
  unit_tiny_matrix_algebra
  tiny_matrix_algebra.cpp
  )

target_link_libraries(
  unit_tiny_matrix_algebra
  tiny
  Boost::unit_test_framework
  )

add_test(
  unit_tiny_matrix_algebra
  unit_tiny_matrix_algebra
  )
