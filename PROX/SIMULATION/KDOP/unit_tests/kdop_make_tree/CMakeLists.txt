add_executable(
  unit_kdop_make_tree
  kdop_make_tree.cpp
  )

target_link_libraries(
  unit_kdop_make_tree
  kdop
  Boost::unit_test_framework
  )

add_test( 
  unit_kdop_make_tree
  unit_kdop_make_tree
  )
