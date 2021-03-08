add_executable(
  unit_util_folder_locator
  util_folder_locator.cpp
  )

target_link_libraries(
  unit_util_folder_locator
  util
  Boost::unit_test_framework
  )

add_test(
  unit_util_folder_locator
  unit_util_folder_locator
  )
