file(GLOB
  GEOMETRY_SOURCES
  ${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/GEOMETRY/GEOMETRY/src/*.cpp
  )

file(GLOB
  GEOMETRY_HEADERS
  ${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/GEOMETRY/GEOMETRY/include/*.h
  )

add_library(geometry
  ${GEOMETRY_SOURCES}
  ${GEOMETRY_HEADERS}
  )

set_target_properties(geometry PROPERTIES PUBLIC_HEADER "${GEOMETRY_HEADERS}")

target_include_directories(geometry
  PUBLIC
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/GEOMETRY/GEOMETRY/include>
  $<INSTALL_INTERFACE:include>
  )

target_link_libraries(geometry
  PUBLIC
  convex
  tiny
  )

install(
  TARGETS geometry
  ARCHIVE DESTINATION lib
  PUBLIC_HEADER DESTINATION include
  )
