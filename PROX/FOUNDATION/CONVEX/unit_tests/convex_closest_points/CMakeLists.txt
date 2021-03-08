add_executable(
  unit_convex_closest_points
  convex_closest_points.cpp
  )

target_link_libraries(
  unit_convex_closest_points
  convex
  Boost::unit_test_framework
  )

add_test(
  unit_convex_closest_points
  unit_convex_closest_points
  )
