add_executable(
  unit_util_timer
  util_timer.cpp
  )

target_link_libraries(
  unit_util_timer
  util
  Boost::unit_test_framework
  )

add_test(
  unit_util_timer
  unit_util_timer
  )
