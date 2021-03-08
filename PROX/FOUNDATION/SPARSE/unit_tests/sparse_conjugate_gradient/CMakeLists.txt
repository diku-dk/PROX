INCLUDE_DIRECTORIES( 
  ${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/SPARSE/SPARSE/include 
  ${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/UTIL/UTIL/include 
  ${Boost_INCLUDE_DIRS}
)

add_executable(
  unit_sparse_conjugate_gradient
  sparse_conjugate_gradient.cpp
  )

target_link_libraries(
  unit_sparse_conjugate_gradient
  util
  Boost::unit_test_framework
  )

add_test(
  unit_sparse_conjugate_gradient
  unit_sparse_conjugate_gradient
  )
