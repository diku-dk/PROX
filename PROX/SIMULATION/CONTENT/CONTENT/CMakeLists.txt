file(GLOB 
  CONTENT_SOURCES 
  src/*.cpp
  )

file(GLOB
  CONTENT_HEADERS
  include/*.h
  )

add_library(content STATIC
  ${CONTENT_SOURCES}
  ${CONTENT_HEADERS}
  )

target_include_directories(content
  PUBLIC
  ${PROJECT_SOURCE_DIR}/PROX/SIMULATION/CONTENT/CONTENT/include
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
 )
