file(GLOB
  MESH_ARRAY_SOURCES
  ${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/MESH_ARRAY/MESH_ARRAY/src/*.cpp
  ${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/MESH_ARRAY/MESH_ARRAY/src/factory/*.cpp
  )
  
file(GLOB
  MESH_ARRAY_HEADERS
  ${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/MESH_ARRAY/MESH_ARRAY/include/*.h
  ${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/MESH_ARRAY/MESH_ARRAY/include/factory/*.h
  )

add_library(mesh_array STATIC
  ${MESH_ARRAY_SOURCES}
  ${MESH_ARRAY_HEADERS}
  )

target_include_directories(mesh_array
  PUBLIC
  ${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/MESH_ARRAY/MESH_ARRAY/include
  PRIVATE
  ${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/MESH_ARRAY/MESH_ARRAY/include/factory
  )

target_link_libraries(mesh_array
  PUBLIC
  tiny
  util
  qhull
  tinyxml
  tetgen
  Boost::boost
  Boost::system
  Boost::filesystem
)

install(
  TARGETS mesh_array
  ARCHIVE DESTINATION lib
 )
