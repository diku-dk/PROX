add_executable(
  unit_tiny_poly_roots
  tiny_polynomial_roots.cpp
  )

target_link_libraries(
  unit_tiny_poly_roots
  tiny
  Boost::unit_test_framework
  )

add_test(
  unit_tiny_poly_roots
  unit_tiny_poly_roots
  )
