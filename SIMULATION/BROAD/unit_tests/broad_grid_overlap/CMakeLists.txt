INCLUDE_DIRECTORIES( 
  ${PROJECT_SOURCE_DIR}/SIMULATION/BROAD/BROAD/include  
  ${Boost_INCLUDE_DIRS}
  )

ADD_EXECUTABLE(
  unit_broad_grid_overlap
  broad_grid_overlap.cpp
  )

TARGET_LINK_LIBRARIES(
  unit_broad_grid_overlap 
  ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
  )

ADD_TEST(
  unit_broad_grid_overlap
  unit_broad_grid_overlap
  )

