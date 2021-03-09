file(GLOB
  UTIL_SOURCES
  ${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/UTIL/UTIL/src/*.cpp 
  )

file(GLOB 
  UTIL_HEADERS 
  ${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/UTIL/UTIL/include/*.h 
  )

add_library(util STATIC
  ${UTIL_SOURCES} 
  ${UTIL_HEADERS} 
  )

target_include_directories(util
  PUBLIC
  ${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/UTIL/UTIL/include 
  )

target_link_libraries(util
  PUBLIC
  Boost::filesystem
  )

install(
  TARGETS util
  ARCHIVE DESTINATION lib
  )
