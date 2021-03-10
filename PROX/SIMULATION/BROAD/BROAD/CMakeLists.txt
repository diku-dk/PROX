file(GLOB
  BROAD_SOURCES
  ${PROJECT_SOURCE_DIR}/PROX/SIMULATION/BROAD/BROAD/src/*.cpp 
  )

file(GLOB
  BROAD_HEADERS
  ${PROJECT_SOURCE_DIR}/PROX/SIMULATION/BROAD/BROAD/include/*.h 
  )

add_library(broad
  ${BROAD_SOURCES} 
  ${BROAD_HEADERS}
  )

set_target_properties(broad PROPERTIES PUBLIC_HEADER "${BROAD_HEADERS}")

target_include_directories(broad
  PUBLIC
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/PROX/SIMULATION/BROAD/BROAD/include>
  $<INSTALL_INTERFACE:include>
  )

target_link_libraries(broad
  PUBLIC
  tiny
  util
  )

install(
  TARGETS broad
  ARCHIVE DESTINATION lib
  PUBLIC_HEADER DESTINATION include
 )

