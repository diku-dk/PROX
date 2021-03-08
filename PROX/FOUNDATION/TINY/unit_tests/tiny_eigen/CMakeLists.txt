add_executable(
  unit_tiny_eigen
  tiny_eigen.cpp
  )

target_link_libraries(
  unit_tiny_eigen
  tiny
  Boost::unit_test_framework
  )

add_test(
  unit_tiny_eigen
  unit_tiny_eigen
  )
