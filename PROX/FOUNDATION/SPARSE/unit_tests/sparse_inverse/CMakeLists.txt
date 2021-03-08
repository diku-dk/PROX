add_executable(
  unit_sparse_inverse
  sparse_inverse.cpp
  )

target_link_libraries(
  unit_sparse_inverse
  sparse
  Boost::unit_test_framework
  )

add_test(
  unit_sparse_inverse
  unit_sparse_inverse
  )
