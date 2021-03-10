file(GLOB
  MESH_ARRAY_SOURCES
  ${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/MESH_ARRAY/MESH_ARRAY/src/*.cpp
  )
  
file(GLOB
  MESH_ARRAY_HEADERS
  ${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/MESH_ARRAY/MESH_ARRAY/include/*.h
  )

add_library(mesh_array
  ${MESH_ARRAY_SOURCES}
  ${MESH_ARRAY_HEADERS}
  )

set_target_properties(mesh_array PROPERTIES PUBLIC_HEADER "${MESH_ARRAY_HEADERS}")

target_include_directories(mesh_array
  PUBLIC
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/MESH_ARRAY/MESH_ARRAY/include>
  $<INSTALL_INTERFACE:include>
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
  PUBLIC_HEADER DESTINATION include
 )
