add_executable(
  unit_prox_solver
  prox_solver.cpp
  )

target_link_libraries(
  unit_prox_solver
  prox
  Boost::unit_test_framework
  )

add_test(
  unit_prox_solver
  unit_prox_solver
  )
