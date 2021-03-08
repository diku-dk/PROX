file(GLOB 
  NARROW_SOURCES
  src/*.cpp
  )

file(GLOB
  NARROW_HEADERS
  include/*.h
  )

add_library(narrow
  ${NARROW_SOURCES}  
  ${NARROW_HEADERS}
  )

target_include_directories(narrow
  PUBLIC
  ${PROJECT_SOURCE_DIR}/PROX/SIMULATION/NARROW/NARROW/include
  )

target_link_libraries(narrow
  PUBLIC
  kdop
  tiny
  geometry
  mesh_array
  convex
  )

install(
  TARGETS narrow
  RUNTIME DESTINATION ${PROJECT_SOURCE_DIR}/lib/}
  COMPONENT library
  )
