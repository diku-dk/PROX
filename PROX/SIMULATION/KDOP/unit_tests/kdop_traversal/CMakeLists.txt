add_executable(
  unit_kdop_traversal
  kdop_traversal.cpp
  )

target_link_libraries(
  unit_kdop_traversal
  kdop
  Boost::unit_test_framework
  )

add_test( 
  unit_kdop_traversal
  unit_kdop_traversal
  )
