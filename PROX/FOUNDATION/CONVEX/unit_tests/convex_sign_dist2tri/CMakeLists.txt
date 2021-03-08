add_executable(
  unit_convex_signed_distance_to_triangle
  convex_sign_dist2tri.cpp
  )

target_link_libraries(
  unit_convex_signed_distance_to_triangle 
  convex
  Boost::unit_test_framework
  )

add_test(
  unit_convex_signed_distance_to_triangle
  unit_convex_signed_distance_to_triangle
  )
