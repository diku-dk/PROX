add_executable(
  unit_convex_motion_interpolation
  convex_motion_interpolation.cpp
  )

target_link_libraries(
  unit_convex_motion_interpolation 
  convex
  Boost::unit_test_framework
  )

add_test(
  unit_convex_motion_interpolation
  unit_convex_motion_interpolation
  )
