file(GLOB
  GL3_SOURCES
  src/*.cpp
  )

file(GLOB
  GL3_HEADERS
  include/*.h
  )

add_library(gl3
  ${GL3_SOURCES}
  ${GL3_HEADERS}
  )

target_include_directories(gl3
  PUBLIC
  ${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/GL3/GL3/include
  ${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/UTIL/UTIL/include
  ${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/TINY/TINY/include
  ${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/GEOMETRY/GEOMETRY/include
  ${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/MESH_ARRAY/MESH_ARRAY/include
  )

target_link_libraries(gl3
  PUBLIC
  tiny
  geometry
  mesh_array
  util
  OpenGL::GL
  GLEW::glew
  glm::glm
  )

install(
  TARGETS gl3
  RUNTIME DESTINATION ${PROJECT_SOURCE_DIR}/lib/}
  COMPONENT library
   )
