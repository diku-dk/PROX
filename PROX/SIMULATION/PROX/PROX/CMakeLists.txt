file(GLOB 
  PROX_SOURCES 
  src/*.cpp
  )

file(GLOB
  PROX_HEADERS
  include/*.h
  include/steppers/*.h
  include/solvers/*.h
  include/solvers/sub/*.h
  include/solvers/strategies/*.h
  )

ADD_LIBRARY( 
  prox 
  ${PROX_SOURCES}  
  ${PROX_HEADERS}
)

target_include_directories(kdop
  PUBLIC
  ${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/TINY/TINY/include 
  ${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/SPARSE/SPARSE/include 
  ${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/GEOMETRY/GEOMETRY/include 
  ${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/MESH_ARRAY/MESH_ARRAY/include
  ${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/CONVEX/CONVEX/include 
  ${PROJECT_SOURCE_DIR}/PROX/SIMULATION/NARROW/NARROW/include
  ${PROJECT_SOURCE_DIR}/PROX/SIMULATION/KDOP/KDOP/include
  ${PROJECT_SOURCE_DIR}/PROX/SIMULATION/BROAD/BROAD/include 
  ${PROJECT_SOURCE_DIR}/PROX/SIMULATION/PROX/PROX/include 
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
  RUNTIME DESTINATION ${PROJECT_SOURCE_DIR}/lib/}
  COMPONENT library
  )


# Add folders to filters

file(GLOB
  PROX_SOLVERS
  include/solvers/*.h
  include/solvers/sub/*.h
  include/solvers/strategies/*.h
  )

source_group(solvers   FILES ${PROX_SOLVERS})

file(GLOB
  PROX_STEPPERS
  include/steppers/*.h
  )

source_group(steppers  FILES ${PROX_STEPPERS})
