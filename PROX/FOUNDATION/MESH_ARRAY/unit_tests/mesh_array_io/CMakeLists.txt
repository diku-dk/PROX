configure_file(
  "${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/MESH_ARRAY/unit_tests/mesh_array_io/mesh_array_data_path.h.in"
  "${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/MESH_ARRAY/unit_tests/mesh_array_io/mesh_array_data_path.h"
  @ONLY
  )

add_executable(
  unit_mesh_array_io 
  mesh_array_io.cpp
  mesh_array_data_path.h
  )

target_link_libraries(
  unit_mesh_array_io
  Boost::unit_test_framework
  mesh_array
  )

add_test( 
  unit_mesh_array_io
  unit_mesh_array_io
  )
