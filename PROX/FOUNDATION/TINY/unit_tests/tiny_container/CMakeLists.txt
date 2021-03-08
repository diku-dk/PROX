INCLUDE_DIRECTORIES( 
  ${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/TINY/TINY/include 
  ${Boost_INCLUDE_DIRS}
  )

add_executable(
  unit_tiny_container
  tiny_container.cpp
  )

target_link_libraries(
  unit_tiny_container
  Boost::unit_test_framework
  )

add_test(
  unit_tiny_container
  unit_tiny_container
  )

