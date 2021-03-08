add_executable(
  unit_sparse_transpose_crm
  sparse_transpose_crm.cpp
  )

target_link_libraries(
  unit_sparse_transpose_crm
  sparse
  Boost::unit_test_framework
  )

add_test(
  unit_sparse_transpose_crm
  unit_sparse_transpose_crm
  )
