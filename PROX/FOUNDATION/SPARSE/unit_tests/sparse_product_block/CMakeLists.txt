INCLUDE_DIRECTORIES( 
  ${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/SPARSE/SPARSE/include 
  ${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/UTIL/UTIL/include 
  ${Boost_INCLUDE_DIRS}
)

add_executable(
  unit_sparse_product_block
  sparse_product_block.cpp
  )

target_link_libraries(
  unit_sparse_product_block
  util
  Boost::unit_test_framework
  )

add_test(
  unit_sparse_product_block
  unit_sparse_product_block
  )
