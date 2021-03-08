INCLUDE_DIRECTORIES( 
  ${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/UTIL/UTIL/include 
  ${Boost_INCLUDE_DIRS}
)

add_executable(
  unit_util_file_locator
  util_file_locator.cpp
  )

target_link_libraries(
  unit_util_file_locator
  ${Boost_SYSTEM_LIBRARY}
  ${Boost_FILESYSTEM_LIBRARY}
  Boost::unit_test_framework
  util
)

add_test(
  unit_util_file_locator
  unit_util_file_locator
  )

