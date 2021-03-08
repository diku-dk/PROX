INCLUDE_DIRECTORIES( 
  ${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/TINY/TINY/include 
  ${Boost_INCLUDE_DIRS}
  )

add_executable(
  unit_tiny_matrix_algebra
  tiny_matrix_algebra.cpp
  )

target_link_libraries(
  unit_tiny_matrix_algebra
  Boost::unit_test_framework
  )

add_test(
  unit_tiny_matrix_algebra
  unit_tiny_matrix_algebra
  )

