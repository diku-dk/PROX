add_executable(
  unit_convex_outside_vertex_edge
  convex_outside_vertex_edge.cpp
  )

target_link_libraries(
  unit_convex_outside_vertex_edge 
  convex
  Boost::unit_test_framework
  )

add_test(
  unit_convex_outside_vertex_edge
  unit_convex_outside_vertex_edge
  )

