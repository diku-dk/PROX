file(GLOB
  MASS_SOURCES
  ${PROJECT_SOURCE_DIR}/PROX/SIMULATION/MASS/MASS/src/*.cpp
  )

file(GLOB
  MASS_HEADERS
  ${PROJECT_SOURCE_DIR}/PROX/SIMULATION/MASS/MASS/include/*.h
  )

add_library(mass STATIC
  ${MASS_SOURCES}
  ${MASS_HEADERS}
  )

target_include_directories(mass
  PUBLIC
  ${PROJECT_SOURCE_DIR}/PROX/SIMULATION/MASS/MASS/include
  )

target_link_libraries(mass
  PUBLIC
  tiny
  )

install(
  TARGETS mass
  ARCHIVE DESTINATION lib
  )
