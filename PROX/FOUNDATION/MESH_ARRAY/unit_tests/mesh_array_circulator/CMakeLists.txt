INCLUDE_DIRECTORIES( 
  ${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/MESH_ARRAY/MESH_ARRAY/include  
  ${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/TINY/TINY/include
  ${Boost_INCLUDE_DIRS}
)

add_executable(
  unit_mesh_array_circulator 
  mesh_array_circulator.cpp
  )

target_link_libraries(
  unit_mesh_array_circulator
  mesh_array
  tiny
  Boost::unit_test_framework
  )

add_test( 
  unit_mesh_array_circulator
  unit_mesh_array_circulator
  )

