add_executable(
  unit_convex_outside_edge_face
  convex_outside_edge_face.cpp
  )

target_link_libraries(
  unit_convex_outside_edge_face 
  convex
  Boost::unit_test_framework
  )

add_test(
  unit_convex_outside_edge_face
  unit_convex_outside_edge_face
  )
