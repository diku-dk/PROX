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

target_include_directories(kdop
  PUBLIC
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

