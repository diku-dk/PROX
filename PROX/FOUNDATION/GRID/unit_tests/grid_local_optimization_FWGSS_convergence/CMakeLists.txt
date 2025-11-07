INCLUDE_DIRECTORIES(
  ${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/GRID/GRID/include
  ${Boost_INCLUDE_DIRS}
  )

ADD_EXECUTABLE(
  unit_grid_local_optimization_FWGSS_convergence
  grid_local_optimization_FWGSS_convergence.cpp
  )

TARGET_LINK_LIBRARIES(
  unit_grid_local_optimization_FWGSS_convergence
  ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
  )

ADD_TEST(
  unit_grid_local_optimization_FWGSS_convergence
    unit_grid_local_optimization_FWGSS_convergence
  )

