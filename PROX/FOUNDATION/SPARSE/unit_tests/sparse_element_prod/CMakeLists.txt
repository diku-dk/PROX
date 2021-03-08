INCLUDE_DIRECTORIES( 
  ${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/SPARSE/SPARSE/include 
  ${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/UTIL/UTIL/include 
  ${Boost_INCLUDE_DIRS}
)

add_executable(
  unit_sparse_element_prod
  sparse_element_prod.cpp
  )

target_link_libraries(
  unit_sparse_element_prod
  util
  Boost::unit_test_framework
  )

add_test(
  unit_sparse_element_prod
  unit_sparse_element_prod
  )
