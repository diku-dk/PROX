INCLUDE_DIRECTORIES(
  ${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/GRID/GRID/include
  ${Boost_INCLUDE_DIRS}
  )

ADD_EXECUTABLE(
  unit_grid_sdf_sdf_ccd_ground_truth_analytical
  grid_sdf_sdf_ccd_ground_truth_analytical.cpp
  )

TARGET_LINK_LIBRARIES(
  unit_grid_sdf_sdf_ccd_ground_truth_analytical
  ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
  )

ADD_TEST(
  unit_grid_sdf_sdf_ccd_ground_truth_analytical
    unit_grid_sdf_sdf_ccd_ground_truth_analytical
  )

