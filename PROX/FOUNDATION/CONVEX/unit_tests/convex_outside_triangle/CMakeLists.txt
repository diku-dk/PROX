add_executable(
  unit_convex_outside_triangle
  convex_outside_tri.cpp
  )

target_link_libraries(
  unit_convex_outside_triangle 
  convex
  Boost::unit_test_framework
  )

add_test(
  unit_convex_outside_triangle
  unit_convex_outside_triangle
  )
