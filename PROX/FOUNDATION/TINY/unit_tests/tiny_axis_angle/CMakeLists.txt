add_executable(
  unit_tiny_axis_angle
  tiny_axis_angle.cpp
  )

target_link_libraries(
  unit_tiny_axis_angle
  tiny
  Boost::unit_test_framework
  )

add_test(
  unit_tiny_axis_angle
  unit_tiny_axis_angle
  )
