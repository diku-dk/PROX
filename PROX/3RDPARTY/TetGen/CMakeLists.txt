
add_library(tetgen STATIC
  predicates.cxx
  tetgen.cxx
  tetgen.h
  )

target_include_directories(tetgen
  PUBLIC
  ${PROJECT_SOURCE_DIR}/PROX/3RDPARTY/TetGen
  )

target_compile_definitions(tetgen
  PUBLIC
  TETLIBRARY
)

install(
  TARGETS tetgen
  ARCHIVE DESTINATION lib
  )


