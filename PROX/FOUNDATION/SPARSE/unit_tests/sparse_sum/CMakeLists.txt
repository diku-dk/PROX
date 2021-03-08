add_executable(
  unit_sparse_sum
  sparse_sum.cpp
  )

target_link_libraries(
  unit_sparse_sum
  sparse
  Boost::unit_test_framework
  )

add_test(
  unit_sparse_sum
  unit_sparse_sum
  )
