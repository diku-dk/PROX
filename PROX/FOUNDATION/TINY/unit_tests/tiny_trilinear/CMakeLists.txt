add_executable(
  unit_tiny_trilinear
  tiny_trilinear.cpp
  )

target_link_libraries(
  unit_tiny_trilinear
  tiny
  Boost::unit_test_framework
  )

add_test(
  unit_tiny_trilinear
  unit_tiny_trilinear
  )
