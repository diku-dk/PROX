INCLUDE_DIRECTORIES( 
  ${PROJECT_SOURCE_DIR}/PROX/SIMULATION/BROAD/BROAD/include  
  ${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/TINY/TINY/include
${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/UTIL/UTIL/include
  ${Boost_INCLUDE_DIRS}
  )

add_executable(
  unit_broad_grid_overlap
  broad_grid_overlap.cpp
  )

target_link_libraries(
  unit_broad_grid_overlap
  tiny
  util
  Boost::unit_test_framework
  )

add_test(
  unit_broad_grid_overlap
  unit_broad_grid_overlap
  )

