INCLUDE_DIRECTORIES(
  ${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/GRID/GRID/include
  ${Boost_INCLUDE_DIRS}
  )

ADD_EXECUTABLE(
  unit_grid_sdf_sdf_ccd_experiments
  grid_sdf_sdf_ccd_experiments.cpp
  )

TARGET_LINK_LIBRARIES(
  unit_grid_sdf_sdf_ccd_experiments
  ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
  )

ADD_TEST(
  unit_grid_sdf_sdf_ccd_experiments
    unit_grid_sdf_sdf_ccd_experiments
  )

