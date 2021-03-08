add_executable(
  unit_narrow_dispatch
  narrow_dispatch.cpp
  )

target_link_libraries(
  unit_narrow_dispatch
  narrow
  Boost::unit_test_framework
  )

add_test(
  unit_narrow_dispatch
  unit_narrow_dispatch
  )
