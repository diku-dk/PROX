add_executable(
  unit_convex_motion_utils
  convex_motion_utils.cpp
  )

target_link_libraries(
  unit_convex_motion_utils
  convex
  Boost::unit_test_framework
  )

add_test(
  unit_convex_motion_utils
  unit_convex_motion_utils
  )
