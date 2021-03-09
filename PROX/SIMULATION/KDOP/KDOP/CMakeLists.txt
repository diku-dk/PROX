file(GLOB 
  KDOP_SOURCES
  src/*.cpp
  )

file(GLOB
  KDOP_HEADERS
  include/*.h
  )

add_library(kdop STATIC
  ${KDOP_SOURCES}  
  ${KDOP_HEADERS}
)

target_include_directories(kdop
  PUBLIC
  ${PROJECT_SOURCE_DIR}/PROX/SIMULATION/KDOP/KDOP/include
  )

target_link_libraries(kdop
  PUBLIC
  tiny
  geometry
  mesh_array
  util
  )

install(
  TARGETS kdop
  ARCHIVE DESTINATION lib
  )
