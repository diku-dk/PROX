add_executable(
  unit_util_file_locator
  util_file_locator.cpp
  )

target_link_libraries(
  unit_util_file_locator
  util
  Boost::unit_test_framework
  )

add_test(
  unit_util_file_locator
  unit_util_file_locator
  )
