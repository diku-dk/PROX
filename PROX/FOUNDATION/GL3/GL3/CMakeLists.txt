file(GLOB
  GL3_SOURCES
  src/*.cpp
  )

file(GLOB
  GL3_HEADERS
  include/*.h
  )

add_library(gl3 STATIC
  ${GL3_SOURCES}
  ${GL3_HEADERS}
  )

target_include_directories(gl3
  PUBLIC
  ${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/GL3/GL3/include
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
  ARCHIVE DESTINATION lib
  )
