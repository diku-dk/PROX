INCLUDE_DIRECTORIES( 
  ${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/MESH_ARRAY/MESH_ARRAY/include  
  ${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/TINY/TINY/include
  ${Boost_INCLUDE_DIRS}
)

add_executable(
  unit_mesh_array_mesh 
  mesh_array_mesh.cpp
  )

target_link_libraries(
  unit_mesh_array_mesh
  Boost::unit_test_framework
  mesh_array
  tiny
  )

add_test( 
  unit_mesh_array_mesh
  unit_mesh_array_mesh
  )

