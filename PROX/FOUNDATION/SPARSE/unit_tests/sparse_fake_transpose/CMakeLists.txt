add_executable(
  unit_sparse_fake_transpose
  sparse_fake_transpose.cpp
  )

target_link_libraries(
  unit_sparse_fake_transpose
  sparse
  Boost::unit_test_framework
  )

add_test(
  unit_sparse_fake_transpose
  unit_sparse_fake_transpose
  )
