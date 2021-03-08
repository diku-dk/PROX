INCLUDE_DIRECTORIES( 
  ${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/TINY/TINY/include 
  ${Boost_INCLUDE_DIRS}
  )

add_executable(
  unit_tiny_precision
  tiny_precision.cpp
  )

target_link_libraries(
  unit_tiny_precision
  Boost::unit_test_framework
  )

add_test(
  unit_tiny_precision
  unit_tiny_precision
  )

