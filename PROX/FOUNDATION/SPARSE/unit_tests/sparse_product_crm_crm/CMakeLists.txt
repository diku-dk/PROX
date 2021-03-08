add_executable(
  unit_sparse_product_crm_crm
  sparse_product_crm_crm.cpp
  )

target_link_libraries(
  unit_sparse_product_crm_crm
  sparse
  Boost::unit_test_framework
  )

add_test(
  unit_sparse_product_crm_crm
  unit_sparse_product_crm_crm
  )
