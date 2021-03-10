file(GLOB
  MASS_SOURCES
  ${PROJECT_SOURCE_DIR}/PROX/SIMULATION/MASS/MASS/src/*.cpp
  )

file(GLOB
  MASS_HEADERS
  ${PROJECT_SOURCE_DIR}/PROX/SIMULATION/MASS/MASS/include/*.h
  )

add_library(mass
  ${MASS_SOURCES}
  ${MASS_HEADERS}
  )

set_target_properties(mass PROPERTIES PUBLIC_HEADER "${MASS_HEADERS}")

target_include_directories(mass
  PUBLIC
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/PROX/SIMULATION/MASS/MASS/include>
  $<INSTALL_INTERFACE:include>
  )

target_link_libraries(mass
  PUBLIC
  tiny
  )

install(
  TARGETS mass
  ARCHIVE DESTINATION lib
  PUBLIC_HEADER DESTINATION include
  )
