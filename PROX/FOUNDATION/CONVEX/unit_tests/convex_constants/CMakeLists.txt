add_executable(
  unit_convex_constants
  convex_constants.cpp
  )

target_link_libraries(
  unit_convex_constants 
  convex
  Boost::unit_test_framework
  )

add_test(
  unit_convex_constants
  unit_convex_constants
  )
