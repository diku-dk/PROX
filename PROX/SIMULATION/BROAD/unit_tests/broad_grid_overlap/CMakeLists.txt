add_executable(
  unit_broad_grid_overlap
  broad_grid_overlap.cpp
  )

target_link_libraries(
  unit_broad_grid_overlap
  broad
  Boost::unit_test_framework
  )

add_test(
  unit_broad_grid_overlap
  unit_broad_grid_overlap
  )
