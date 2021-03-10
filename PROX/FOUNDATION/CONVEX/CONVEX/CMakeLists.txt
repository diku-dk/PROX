file(GLOB
  CONVEX_SOURCES
  ${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/CONVEX/CONVEX/src/*.cpp 
  )

file(GLOB 
  CONVEX_HEADERS 
  ${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/CONVEX/CONVEX/include/*.h 
  )

add_library(convex
  ${CONVEX_SOURCES} 
  ${CONVEX_HEADERS} 
  )

set_target_properties(convex PROPERTIES PUBLIC_HEADER "${CONVEX_HEADERS}")

target_include_directories(convex
  PUBLIC
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/CONVEX/CONVEX/include>
  $<INSTALL_INTERFACE:include>
  )

target_link_libraries( convex
  PUBLIC
  tiny
  geometry
  )

install(
  TARGETS convex
  ARCHIVE DESTINATION lib
  PUBLIC_HEADER DESTINATION include
  )
