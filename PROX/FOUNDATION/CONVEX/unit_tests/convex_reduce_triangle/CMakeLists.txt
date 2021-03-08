add_executable(
  unit_convex_reduce_triangle
  convex_reduce_triangle.cpp
  )

target_link_libraries(
  unit_convex_reduce_triangle 
  convex
  Boost::unit_test_framework
  )

add_test(
  unit_convex_reduce_triangle
  unit_convex_reduce_triangle
  )
