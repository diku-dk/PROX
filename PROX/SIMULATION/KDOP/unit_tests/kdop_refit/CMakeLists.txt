add_executable(
  unit_kdop_refit
  kdop_refit.cpp
  )

target_link_libraries(
  unit_kdop_refit
  kdop
  Boost::unit_test_framework
  )

add_test( 
  unit_kdop_refit
  unit_kdop_refit
  )
