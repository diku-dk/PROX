add_executable(
  unit_convex_reduce_tetrahedron
  convex_reduce_tetrahedron.cpp
  )

target_link_libraries(
  unit_convex_reduce_tetrahedron 
  convex
  Boost::unit_test_framework
  )

add_test(
  unit_convex_reduce_tetrahedron
  unit_convex_reduce_tetrahedron
  )
