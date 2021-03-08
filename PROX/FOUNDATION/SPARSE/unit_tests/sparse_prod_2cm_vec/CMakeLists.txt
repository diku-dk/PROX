add_executable(
  unit_sparse_prod_2cm_vec
  sparse_prod_2cm_vec.cpp
  )

target_link_libraries(
  unit_sparse_prod_2cm_vec
  sparse
  Boost::unit_test_framework
  )

add_test(
  unit_sparse_prod_2cm_vec
  unit_sparse_prod_2cm_vec
  )
