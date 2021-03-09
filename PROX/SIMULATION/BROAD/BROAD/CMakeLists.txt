file(GLOB
  BROAD_SOURCES
  ${PROJECT_SOURCE_DIR}/PROX/SIMULATION/BROAD/BROAD/src/*.cpp 
  )

file(GLOB
  BROAD_HEADERS
  ${PROJECT_SOURCE_DIR}/PROX/SIMULATION/BROAD/BROAD/include/*.h 
  )

add_library(broad STATIC
  ${BROAD_SOURCES} 
  ${BROAD_HEADERS}
  )

target_include_directories(broad
  PUBLIC
  ${PROJECT_SOURCE_DIR}/PROX/SIMULATION/BROAD/BROAD/include
  )

target_link_libraries(broad
  PUBLIC
  tiny
  util
  )

install(
  TARGETS broad
  ARCHIVE DESTINATION lib
 )

