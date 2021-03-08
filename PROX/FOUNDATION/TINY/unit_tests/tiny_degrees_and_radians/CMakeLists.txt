INCLUDE_DIRECTORIES( 
  ${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/TINY/TINY/include 
  ${Boost_INCLUDE_DIRS}
  )

add_executable(
  unit_tiny_degrees_and_radians
  tiny_degrees_and_radians.cpp
  )

target_link_libraries(
  unit_tiny_degrees_and_radians
  Boost::unit_test_framework
  )

add_test(
  unit_tiny_degrees_and_radians
  unit_tiny_degrees_and_radians
  )
