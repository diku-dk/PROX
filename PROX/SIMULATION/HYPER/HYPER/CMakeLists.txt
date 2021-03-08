file(GLOB
  HYPER_SOURCES
  ${PROJECT_SOURCE_DIR}/PROX/SIMULATION/HYPER/HYPER/src/*.cpp
  )

file(GLOB
  HYPER_HEADERS ${PROJECT_SOURCE_DIR}/PROX/SIMULATION/HYPER/HYPER/include/*.h
  )

add_library(hyper
  ${HYPER_SOURCES}
  ${HYPER_HEADERS}
  )

target_include_directories(hyper
  PUBLIC
  ${PROJECT_SOURCE_DIR}/PROX/SIMULATION/HYPER/HYPER/include
  )

target_link_libraries(hyper
  PUBLIC
  tiny
  sparse
  mesh_array
  kdop
  )

install(
  TARGETS hyper
  RUNTIME DESTINATION ${PROJECT_SOURCE_DIR}/lib/}
  COMPONENT library
 )
