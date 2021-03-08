add_executable(
  unit_convex_simplex_reduce
  convex_simplex_reduce.cpp
  )

target_link_libraries(
  unit_convex_simplex_reduce
  convex
  Boost::unit_test_framework
  )

add_test(
  unit_convex_simplex_reduce
  unit_convex_simplex_reduce
  )
