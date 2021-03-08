add_executable(
  unit_convex_shapes
  convex_shapes.cpp
  )

target_link_libraries(
  unit_convex_shapes 
  convex
  Boost::unit_test_framework
  )

add_test(
  unit_convex_shapes
  unit_convex_shapes
  )
