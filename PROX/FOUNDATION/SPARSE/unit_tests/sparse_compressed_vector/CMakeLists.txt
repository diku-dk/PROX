add_executable(
  unit_sparse_compressed_vector
  sparse_compressed_vector.cpp
  )

target_link_libraries(
  unit_sparse_compressed_vector
  sparse
  Boost::unit_test_framework
  )

add_test(
  unit_sparse_compressed_vector
  unit_sparse_compressed_vector
  )
