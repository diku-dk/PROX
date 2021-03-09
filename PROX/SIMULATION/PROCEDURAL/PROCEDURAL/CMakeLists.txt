file(GLOB 
  PROCEDURAL_SOURCES
  src/*.cpp
  )

file(GLOB
  PROCEDURAL_HEADERS
  include/*.h
  )

add_library(procedural STATIC
  ${PROCEDURAL_SOURCES}  
  ${PROCEDURAL_HEADERS}
)

target_include_directories(procedural
  PUBLIC
  ${PROJECT_SOURCE_DIR}/PROX/SIMULATION/PROCEDURAL/PROCEDURAL/include
  )

target_link_libraries(procedural
  PUBLIC
  util
  mass
  tiny
  mesh_array
  content
  )

install(
  TARGETS procedural
  ARCHIVE DESTINATION lib
  )
