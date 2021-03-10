file(GLOB 
  KDOP_SOURCES
  src/*.cpp
  )

file(GLOB
  KDOP_HEADERS
  include/*.h
  )

add_library(kdop
  ${KDOP_SOURCES}  
  ${KDOP_HEADERS}
)

set_target_properties(kdop PROPERTIES PUBLIC_HEADER "${KDOP_HEADERS}")

target_include_directories(kdop
  PUBLIC
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/PROX/SIMULATION/KDOP/KDOP/include>
  $<INSTALL_INTERFACE:include>
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
  PUBLIC_HEADER DESTINATION include
  )
