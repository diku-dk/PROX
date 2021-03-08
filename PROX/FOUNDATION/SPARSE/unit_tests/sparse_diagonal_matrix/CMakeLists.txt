add_executable(
  unit_sparse_diagonal_matrix
  sparse_diagonal_matrix.cpp
  )

target_link_libraries(
  unit_sparse_diagonal_matrix
  sparse
  Boost::unit_test_framework
  )

add_test(
  unit_sparse_diagonal_matrix
  unit_sparse_diagonal_matrix
  )
