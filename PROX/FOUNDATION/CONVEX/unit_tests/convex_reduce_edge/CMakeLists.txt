add_executable(
  unit_convex_reduce_edge
  convex_reduce_edge.cpp
  )

target_link_libraries(
  unit_convex_reduce_edge 
  convex
  Boost::unit_test_framework
  )

add_test(
  unit_convex_reduce_edge
  unit_convex_reduce_edge
  )
