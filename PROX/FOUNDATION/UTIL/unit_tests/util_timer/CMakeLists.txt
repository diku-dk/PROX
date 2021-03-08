INCLUDE_DIRECTORIES( 
  ${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/UTIL/UTIL/include 
  ${Boost_INCLUDE_DIRS}
)

add_executable(
  unit_util_timer
  util_timer.cpp
  )

target_link_libraries(
  unit_util_timer
  Boost::unit_test_framework
  util
)

add_test(
  unit_util_timer
  unit_util_timer
  )

