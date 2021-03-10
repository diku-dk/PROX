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

set_target_properties(narrow PROPERTIES PUBLIC_HEADER "${NARROW_HEADERS}")

target_include_directories(narrow
  PUBLIC
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/PROX/SIMULATION/NARROW/NARROW/include>
  $<INSTALL_INTERFACE:include>
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
  ARCHIVE DESTINATION lib
  PUBLIC_HEADER DESTINATION include
  )
