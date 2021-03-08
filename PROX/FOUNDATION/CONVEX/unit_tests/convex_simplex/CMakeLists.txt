add_executable(
  unit_convex_simplex
  convex_simplex.cpp
  )

target_link_libraries(
  unit_convex_simplex
  convex
  Boost::unit_test_framework
  )

add_test(
  unit_convex_simplex
  unit_convex_simplex
  )
