INCLUDE_DIRECTORIES(
  ${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/GRID/GRID/include
  ${Boost_INCLUDE_DIRS}
  )

ADD_EXECUTABLE(
  unit_grid_test_grid
  grid_test_grid.cpp
  )

TARGET_LINK_LIBRARIES(
  unit_grid_test_grid
  ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
  )

ADD_TEST(
  unit_grid_test_grid
  unit_grid_test_grid
  )

