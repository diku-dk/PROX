configure_file(
  "${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/MESH_ARRAY/unit_tests/mesh_array_tetgen/mesh_array_data_path.h.in"
  "${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/MESH_ARRAY/unit_tests/mesh_array_tetgen/mesh_array_data_path.h"
  @ONLY
  )

add_executable(
  unit_mesh_array_tetgen 
  mesh_array_tetgen.cpp
  mesh_array_data_path.h
  )

target_link_libraries(
  unit_mesh_array_tetgen
  mesh_array
  Boost::unit_test_framework
  )

add_test( 
  unit_mesh_array_tetgen
  unit_mesh_array_tetgen
  )
