add_executable(
  unit_tiny_coordsys_prod
  tiny_coordsys_prod.cpp
  )

target_link_libraries(
  unit_tiny_coordsys_prod
  tiny
  Boost::unit_test_framework
  )

add_test(
  unit_tiny_coordsys_prod
  unit_tiny_coordsys_prod
  )
