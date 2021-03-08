add_executable(
  unit_sparse_row_product
  sparse_row_product.cpp
  )

target_link_libraries(
  unit_sparse_row_product
  sparse
  Boost::unit_test_framework
  )

add_test(
  unit_sparse_row_product
  unit_sparse_row_product
  )
