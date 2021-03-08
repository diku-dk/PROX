add_executable(
  unit_tiny_quaternion_rotate
  tiny_quaternion_rotate.cpp
  )

target_link_libraries(
  unit_tiny_quaternion_rotate
  tiny
  Boost::unit_test_framework
  )

add_test(
  unit_tiny_quaternion_rotate
  unit_tiny_quaternion_rotate
  )
