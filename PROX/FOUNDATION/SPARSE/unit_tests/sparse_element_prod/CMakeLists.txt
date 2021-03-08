add_executable(
  unit_sparse_element_prod
  sparse_element_prod.cpp
  )

target_link_libraries(
  unit_sparse_element_prod
  sparse
  Boost::unit_test_framework
  )

add_test(
  unit_sparse_element_prod
  unit_sparse_element_prod
  )
