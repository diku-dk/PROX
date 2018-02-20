INCLUDE_DIRECTORIES( 
  ${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/UTIL/UTIL/include 
  ${Boost_INCLUDE_DIRS}
)

ADD_EXECUTABLE(
  unit_util_folder_locator
  util_folder_locator.cpp
  )

TARGET_LINK_LIBRARIES(
  unit_util_folder_locator
  ${Boost_SYSTEM_LIBRARY}
  ${Boost_FILESYSTEM_LIBRARY}
  ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
  util
)

ADD_TEST(
  unit_util_folder_locator
  unit_util_folder_locator
  )

