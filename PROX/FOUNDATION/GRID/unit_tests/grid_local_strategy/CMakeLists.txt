INCLUDE_DIRECTORIES(
  ${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/GRID/GRID/include
  ${Boost_INCLUDE_DIRS}
  )

ADD_EXECUTABLE(
  unit_grid_local_strategy
  grid_local_strategy.cpp
  )

TARGET_LINK_LIBRARIES(
  unit_grid_local_strategy
  ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
  )

ADD_TEST(
  unit_grid_local_strategy
  unit_grid_local_strategy
  )

