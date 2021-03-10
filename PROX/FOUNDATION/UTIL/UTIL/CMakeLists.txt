file(GLOB
  UTIL_SOURCES
  ${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/UTIL/UTIL/src/*.cpp 
  )

file(GLOB 
  UTIL_HEADERS 
  ${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/UTIL/UTIL/include/*.h 
  )

add_library(util
  ${UTIL_SOURCES} 
  ${UTIL_HEADERS} 
  )

set_target_properties(util PROPERTIES PUBLIC_HEADER "${UTIL_HEADERS}")

target_include_directories(util
  PUBLIC
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/UTIL/UTIL/include>
  $<INSTALL_INTERFACE:include>
  )

target_link_libraries(util
  PUBLIC
  Boost::filesystem
  )

install(
  TARGETS util
  ARCHIVE DESTINATION lib
  PUBLIC_HEADER DESTINATION include
  )
