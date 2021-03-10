file(GLOB 
  PROX_SOURCES 
  src/*.cpp
  )

file(GLOB
  PROX_HEADERS
  include/*.h
  )

add_library(prox
  ${PROX_SOURCES}  
  ${PROX_HEADERS}
)

set_target_properties(prox PROPERTIES PUBLIC_HEADER "${PROX_HEADERS}")

target_include_directories(kdop
  PUBLIC
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/PROX/SIMULATION/PROX/PROX/include>
  $<INSTALL_INTERFACE:include>
  )

target_link_libraries(prox
  PUBLIC
  tiny
  sparse
  geometry
  convex
  mesh_array
  broad
  narrow
  kdop
  )

install(
  TARGETS prox
  ARCHIVE DESTINATION lib
  PUBLIC_HEADER DESTINATION include  
  )
