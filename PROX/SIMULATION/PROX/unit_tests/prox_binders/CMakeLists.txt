add_executable(
  unit_prox_binders
  prox_binders.cpp
  )

target_link_libraries(
  unit_prox_binders
  prox
  Boost::unit_test_framework
  )

add_test(
  unit_prox_binders
  unit_prox_binders
  )
