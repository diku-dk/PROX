add_executable(
  unit_convex_signed_distance_to_vertex_edge_vp
  convex_sign_dist2vert_edge_vp.cpp
  )

target_link_libraries(
  unit_convex_signed_distance_to_vertex_edge_vp 
  convex
  Boost::unit_test_framework
  )

add_test(
  unit_convex_signed_distance_to_vertex_edge_vp
  unit_convex_signed_distance_to_vertex_edge_vp
  )
