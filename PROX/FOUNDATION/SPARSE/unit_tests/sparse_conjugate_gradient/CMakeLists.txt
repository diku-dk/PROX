add_executable(
  unit_sparse_conjugate_gradient
  sparse_conjugate_gradient.cpp
  )

target_link_libraries(
  unit_sparse_conjugate_gradient
  sparse
  Boost::unit_test_framework
  )

add_test(
  unit_sparse_conjugate_gradient
  unit_sparse_conjugate_gradient
  )
