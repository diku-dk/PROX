INCLUDE_DIRECTORIES(
  ${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/GRID/GRID/include
  ${Boost_INCLUDE_DIRS}
  )

ADD_EXECUTABLE(
  unit_grid_local_optimization_FWGSS
  grid_local_optimization_FWGSS.cpp
  )

TARGET_LINK_LIBRARIES(
  unit_grid_local_optimization_FWGSS
  ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
  )

ADD_TEST(
  unit_grid_local_optimization_FWGSS
    unit_grid_local_optimization_FWGSS
  )

