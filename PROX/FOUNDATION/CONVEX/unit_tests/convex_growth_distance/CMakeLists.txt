add_executable(
  unit_convex_growth_distance
  convex_growth_distance.cpp
  )

target_link_libraries(
  unit_convex_growth_distance 
  convex
  Boost::unit_test_framework
  )

add_test(
  unit_convex_growth_distance
  unit_convex_growth_distance
  )
