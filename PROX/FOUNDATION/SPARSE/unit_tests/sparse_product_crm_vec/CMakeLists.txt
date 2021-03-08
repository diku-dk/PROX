add_executable(
  unit_sparse_product_crm_vec
  sparse_product_crm_vec.cpp
  )

target_link_libraries(
  unit_sparse_product_crm_vec
  sparse
  Boost::unit_test_framework
  )

add_test(
  unit_sparse_product_crm_vec
  unit_sparse_product_crm_vec
  )
