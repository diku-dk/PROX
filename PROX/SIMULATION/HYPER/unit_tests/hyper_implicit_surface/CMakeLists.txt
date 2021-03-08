add_executable(
  unit_hyper_implicit_surface
  hyper_implicit_surface.cpp
  )

target_link_libraries(
  unit_hyper_implicit_surface
  hyper
  Boost::unit_test_framework
  )

add_test(
  unit_hyper_implicit_surface
  unit_hyper_implicit_surface
  )
