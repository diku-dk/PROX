add_executable(
  unit_tiny_vector_algebra
  tiny_vector_algebra.cpp
  )

target_link_libraries(
  unit_tiny_vector_algebra
  tiny
  Boost::unit_test_framework
  )

add_test(
  unit_tiny_vector_algebra
  unit_tiny_vector_algebra
  )
