INCLUDE_DIRECTORIES( 
  ${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/SPARSE/SPARSE/include 
  ${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/UTIL/UTIL/include 
  ${Boost_INCLUDE_DIRS}
)

add_executable(
  unit_sparse_diagonal_matrix
  sparse_diagonal_matrix.cpp
  )

target_link_libraries(
  unit_sparse_diagonal_matrix
  util
  Boost::unit_test_framework
  )

add_test(
  unit_sparse_diagonal_matrix
  unit_sparse_diagonal_matrix
  )
