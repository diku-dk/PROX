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

set_target_properties(hyper PROPERTIES PUBLIC_HEADER "${HYPER_HEADERS}")

target_include_directories(hyper
  PUBLIC
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/PROX/SIMULATION/HYPER/HYPER/include>
  $<INSTALL_INTERFACE:include>
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
  ARCHIVE DESTINATION lib
  PUBLIC_HEADER DESTINATION include
)
