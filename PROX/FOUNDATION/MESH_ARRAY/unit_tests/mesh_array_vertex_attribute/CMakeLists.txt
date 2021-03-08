INCLUDE_DIRECTORIES(
  ${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/MESH_ARRAY/MESH_ARRAY/include
  ${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/TINY/TINY/include
  ${Boost_INCLUDE_DIRS}
)

add_executable(
  unit_mesh_array_vertex_attribute
  mesh_array_vertex_attribute.cpp
  )

target_link_libraries(
  unit_mesh_array_vertex_attribute
  mesh_array
  tiny
  Boost::unit_test_framework
  )

add_test(
  unit_mesh_array_vertex_attribute
  unit_mesh_array_vertex_attribute
  )

