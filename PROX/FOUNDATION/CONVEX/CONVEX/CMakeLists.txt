file(GLOB
  CONVEX_SOURCES
  ${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/CONVEX/CONVEX/src/*.cpp 
  )

file(GLOB 
  CONVEX_HEADERS 
  ${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/CONVEX/CONVEX/include/*.h 
  )

add_library(convex STATIC
  ${CONVEX_SOURCES} 
  ${CONVEX_HEADERS} 
  )

target_include_directories(convex
  PUBLIC
  ${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/CONVEX/CONVEX/include    
  )

target_link_libraries( convex
  PUBLIC
  tiny
  geometry
  )

install(
  TARGETS convex
  ARCHIVE DESTINATION lib
  )
