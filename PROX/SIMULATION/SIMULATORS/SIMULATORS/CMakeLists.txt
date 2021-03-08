file(GLOB 
  SIMULATORS_SOURCES
  src/*.cpp
  )

file(GLOB
  SIMULATORS_HEADERS
  include/*.h
  private/*.h
  )

add_library(simulators 
  ${SIMULATORS_SOURCES}  
  ${SIMULATORS_HEADERS}
  )


target_include_directories(simulators
  PUBLIC
  ${PROJECT_SOURCE_DIR}/PROX/SIMULATION/SIMULATORS/SIMULATORS/include
  PRIVATE
  ${PROJECT_SOURCE_DIR}/PROX/SIMULATION/SIMULATORS/SIMULATORS/private
  )

target_link_libraries(simulators
  PUBLIC
  content
  tiny
  sparse
  convex
  geometry
  mesh_array
  kdop
  narrow
  broad
  prox
  mass
  util
  tetgen
  tinyxml
  qhull
  Boost::filesystem
  )

install(
  TARGETS simulators
  RUNTIME DESTINATION ${PROJECT_SOURCE_DIR}/lib/}
  COMPONENT library
  )
