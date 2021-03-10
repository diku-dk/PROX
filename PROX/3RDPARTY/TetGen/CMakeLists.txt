
add_library(tetgen
  predicates.cxx
  tetgen.cxx
  tetgen.h
  )

set_target_properties(tetgen PROPERTIES PUBLIC_HEADER "tetgen.h")

target_include_directories(tetgen
  PUBLIC
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/PROX/3RDPARTY/TetGen>
  $<INSTALL_INTERFACE:include>
  )

target_compile_definitions(tetgen
  PUBLIC
  TETLIBRARY
)

install(
  TARGETS tetgen
  ARCHIVE DESTINATION lib
  PUBLIC_HEADER DESTINATION include
  )


