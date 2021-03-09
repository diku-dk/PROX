file(GLOB
  SPARSE_SOURCES
  ${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/SPARSE/SPARSE/src/*.cpp
  )

file(GLOB
  SPARSE_HEADERS
  ${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/SPARSE/SPARSE/include/*.h 
  ${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/SPARSE/SPARSE/include/io/*.h 
  )

add_library(sparse STATIC
  ${SPARSE_SOURCES}
  ${SPARSE_HEADERS}
  )

target_include_directories(sparse
  PUBLIC
  ${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/SPARSE/SPARSE/include
  )

target_link_libraries(sparse
  PUBLIC
  util
  )

install(
  TARGETS sparse
  ARCHIVE DESTINATION lib
 )
