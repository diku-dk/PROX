add_executable(
  unit_tiny_container
  tiny_container.cpp
  )

target_link_libraries(
  unit_tiny_container
  tiny
  Boost::unit_test_framework
  )

add_test(
  unit_tiny_container
  unit_tiny_container
  )
