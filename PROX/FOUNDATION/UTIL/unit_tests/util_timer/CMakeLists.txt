INCLUDE_DIRECTORIES( 
  ${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/UTIL/UTIL/include 
  ${Boost_INCLUDE_DIRS}
)

ADD_EXECUTABLE(
  unit_util_timer
  util_timer.cpp
  )

TARGET_LINK_LIBRARIES(
  unit_util_timer
  ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
  util
)

ADD_TEST(
  unit_util_timer
  unit_util_timer
  )

