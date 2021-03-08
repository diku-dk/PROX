add_executable(
  unit_sparse_transpose_block
  sparse_transpose_block.cpp
  )

target_link_libraries(
  unit_sparse_transpose_block
  sparse
  Boost::unit_test_framework
  )

add_test(
  unit_sparse_transpose_block
  unit_sparse_transpose_block
  )
