add_executable(
  unit_kdop_raycast
  kdop_raycast.cpp
  )

target_link_libraries(
  unit_kdop_raycast
  kdop
  Boost::unit_test_framework
  )

add_test( 
  unit_kdop_raycast
  unit_kdop_raycast
  )
