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

set_target_properties(gl3 PROPERTIES PUBLIC_HEADER "${GL3_HEADERS}")

target_include_directories(gl3
  PUBLIC
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/GL3/GL3/include>
  $<INSTALL_INTERFACE:include>
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
  PUBLIC_HEADER DESTINATION include
  )
