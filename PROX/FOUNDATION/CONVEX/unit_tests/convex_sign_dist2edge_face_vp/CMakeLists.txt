add_executable(
  unit_convex_signed_distance_to_edge_face_vp
  convex_sign_dist2edge_face_vp.cpp
  )

target_link_libraries(
  unit_convex_signed_distance_to_edge_face_vp 
  convex
  Boost::unit_test_framework
  )

add_test(
  unit_convex_signed_distance_to_edge_face_vp
  unit_convex_signed_distance_to_edge_face_vp
  )
