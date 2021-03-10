file(GLOB
  SPARSE_SOURCES
  ${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/SPARSE/SPARSE/src/*.cpp
  )

file(GLOB
  SPARSE_HEADERS
  ${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/SPARSE/SPARSE/include/*.h 
  )

add_library(sparse
  ${SPARSE_SOURCES}
  ${SPARSE_HEADERS}
  )

set_target_properties(sparse PROPERTIES PUBLIC_HEADER "${SPARSE_HEADERS}")

target_include_directories(sparse
  PUBLIC
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/SPARSE/SPARSE/include>
  $<INSTALL_INTERFACE:include>
  )

target_link_libraries(sparse
  PUBLIC
  util
  )

install(
  TARGETS sparse
  ARCHIVE DESTINATION lib
  PUBLIC_HEADER DESTINATION include
 )
