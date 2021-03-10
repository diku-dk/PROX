file(GLOB
  TINY_SOURCES
  ${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/TINY/TINY/src/*.cpp
  )

file(GLOB
  TINY_HEADERS
  ${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/TINY/TINY/include/*.h
  )

add_library(tiny
  ${TINY_SOURCES}
  ${TINY_HEADERS}
  )

set_target_properties(tiny PROPERTIES PUBLIC_HEADER "${TINY_HEADERS}")

target_include_directories(tiny
  PUBLIC
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/TINY/TINY/include>
  $<INSTALL_INTERFACE:include>
  )

target_link_libraries(tiny
  PUBLIC
  Boost::boost
  )

install(
  TARGETS tiny
  ARCHIVE DESTINATION lib
  PUBLIC_HEADER DESTINATION include
  )
 
