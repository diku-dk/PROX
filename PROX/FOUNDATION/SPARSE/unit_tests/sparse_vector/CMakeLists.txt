add_executable(
  unit_sparse_vector
  sparse_vector.cpp
  )

target_link_libraries(
  unit_sparse_vector
  sparse
  Boost::unit_test_framework
  )

add_test(
  unit_sparse_vector
  unit_sparse_vector
  )
