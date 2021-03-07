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

target_include_directories(content
  PUBLIC
  ${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/TINY/TINY/include
  ${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/UTIL/UTIL/include
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
  RUNTIME DESTINATION ${PROJECT_SOURCE_DIR}/lib/}
  COMPONENT library
 )
