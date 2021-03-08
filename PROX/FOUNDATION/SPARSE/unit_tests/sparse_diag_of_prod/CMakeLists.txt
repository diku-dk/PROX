add_executable(
  unit_sparse_diag_of_prod
  sparse_diag_of_prod.cpp
  )

target_link_libraries(
  unit_sparse_diag_of_prod
  sparse
  Boost::unit_test_framework
  )

add_test(
  unit_sparse_diag_of_prod
  unit_sparse_diag_of_prod
  )
