INCLUDE_DIRECTORIES(
  ${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/GRID/GRID/include
  ${Boost_INCLUDE_DIRS}
  )

ADD_EXECUTABLE(
  unit_grid_local_optimization_test_GSS_minimize
  grid_local_optimization_test_GSS_minimize.cpp
  )

TARGET_LINK_LIBRARIES(
  unit_grid_local_optimization_test_GSS_minimize
  ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
  )

ADD_TEST(
  unit_grid_local_optimization_test_GSS_minimize
    unit_grid_local_optimization_test_GSS_minimize
  )

