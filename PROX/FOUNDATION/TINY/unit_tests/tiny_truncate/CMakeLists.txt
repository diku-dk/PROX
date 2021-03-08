add_executable(
  unit_tiny_truncate
  tiny_truncate.cpp
  )

target_link_libraries(
  unit_tiny_truncate
  tiny
  Boost::unit_test_framework
  )

add_test(
  unit_tiny_truncate
  unit_tiny_truncate
  )
