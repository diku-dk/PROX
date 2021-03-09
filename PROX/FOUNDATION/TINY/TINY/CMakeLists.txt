file(GLOB
  TINY_SOURCES
  ${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/TINY/TINY/src/*.cpp
  )

file(GLOB
  TINY_HEADERS
  ${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/TINY/TINY/include/*.h
  )

add_library(tiny STATIC
  ${TINY_SOURCES}
  ${TINY_HEADERS}
  )

target_include_directories(tiny
  PUBLIC
  ${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/TINY/TINY/include
  )

target_link_libraries(tiny
  PUBLIC
  Boost::boost
  )

install(
  TARGETS tiny
  ARCHIVE DESTINATION lib
  )
 
