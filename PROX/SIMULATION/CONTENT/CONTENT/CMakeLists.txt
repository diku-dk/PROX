file(GLOB 
  CONTENT_SOURCES 
  src/*.cpp
  )

file(GLOB
  CONTENT_HEADERS
  include/*.h
  )

add_library(content
  ${CONTENT_SOURCES}
  ${CONTENT_HEADERS}
  )

set_target_properties(content PROPERTIES PUBLIC_HEADER "${CONTENT_HEADERS}")

target_include_directories(content
  PUBLIC
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/PROX/SIMULATION/CONTENT/CONTENT/include>
  $<INSTALL_INTERFACE:include>
 )

target_link_libraries(content
  PUBLIC
  tinyxml
  util
  tiny
  )

install(
  TARGETS content
  ARCHIVE DESTINATION lib
  PUBLIC_HEADER DESTINATION include
 )
