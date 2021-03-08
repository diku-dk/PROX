add_executable(
  unit_sparse_two_column_matrix
  sparse_two_column_matrix.cpp
  )

target_link_libraries(
  unit_sparse_two_column_matrix
  sparse
  Boost::unit_test_framework
  )

add_test(
  unit_sparse_two_column_matrix
  unit_sparse_two_column_matrix
  )
