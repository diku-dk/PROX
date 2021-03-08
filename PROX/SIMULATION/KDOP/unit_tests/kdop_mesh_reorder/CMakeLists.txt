add_executable(
  unit_kdop_mesh_reorder
  kdop_mesh_reorder.cpp
  )

target_link_libraries(
  unit_kdop_mesh_reorder
  kdop
  Boost::unit_test_framework
  )

add_test( 
  unit_kdop_mesh_reorder
  unit_kdop_mesh_reorder
  )
