add_executable(
  unit_sparse_io
  sparse_io.cpp
  )

target_link_libraries(
  unit_sparse_io
  sparse
  Boost::unit_test_framework
  )

add_test(
  unit_sparse_io
  unit_sparse_io
  )
