add_executable(
  unit_sparse_block
  sparse_block.cpp
  )

target_link_libraries(
  unit_sparse_block
  sparse
  Boost::unit_test_framework
  )

add_test(
  unit_sparse_block
  unit_sparse_block
  )
