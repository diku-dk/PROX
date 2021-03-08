add_executable(
  unit_convex_get_scale
  convex_get_scale.cpp
  )

target_link_libraries(
  unit_convex_get_scale 
  convex
  Boost::unit_test_framework
  )

add_test(
  unit_convex_get_scale
  unit_convex_get_scale
  )
