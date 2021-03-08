CONFIGURE_FILE(
  "${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/MESH_ARRAY/unit_tests/mesh_array_tetgen/mesh_array_data_path.h.in"
  "${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/MESH_ARRAY/unit_tests/mesh_array_tetgen/mesh_array_data_path.h"
  @ONLY
  )

INCLUDE_DIRECTORIES( 
  ${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/MESH_ARRAY/MESH_ARRAY/include  
  ${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/TINY/TINY/include
  ${Boost_INCLUDE_DIRS}
)

add_executable(
  unit_mesh_array_tetgen 
  mesh_array_tetgen.cpp
  mesh_array_data_path.h
  )

target_link_libraries(
  unit_mesh_array_tetgen
  Boost::unit_test_framework
  mesh_array
  tiny
  )

add_test( 
  unit_mesh_array_tetgen
  unit_mesh_array_tetgen
  )

