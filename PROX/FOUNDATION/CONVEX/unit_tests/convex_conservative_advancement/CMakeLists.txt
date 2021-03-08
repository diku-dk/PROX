add_executable(
  unit_convex_conservative_advancement
  convex_conservative_advancement.cpp
  )

target_link_libraries(
  unit_convex_conservative_advancement
  convex
  Boost::unit_test_framework
  )

add_test(
  unit_convex_conservative_advancement
  unit_convex_conservative_advancement
  )
