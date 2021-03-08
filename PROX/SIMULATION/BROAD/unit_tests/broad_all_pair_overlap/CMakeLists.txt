add_executable(
  unit_broad_all_pair_overlap
  broad_all_pair_overlap.cpp
  )

target_link_libraries(
  unit_broad_all_pair_overlap
  broad
  Boost::unit_test_framework
  )

add_test(
  unit_broad_all_pair_overlap
  unit_broad_all_pair_overlap
  )
