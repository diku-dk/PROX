add_executable(
  unit_tiny_euler_angles
  tiny_euler_angles.cpp
  )

target_link_libraries(
  unit_tiny_euler_angles
  tiny
  Boost::unit_test_framework
  )

add_test(
  unit_tiny_euler_angles
  unit_tiny_euler_angles
  )
