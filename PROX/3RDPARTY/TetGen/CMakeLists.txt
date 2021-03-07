
add_library(tetgen
  predicates.cxx
  tetgen.cxx
  tetgen.h
  )

target_compile_definitions(tetgen
  PUBLIC
  TETLIBRARY
)

install(
  TARGETS tetgen
  RUNTIME DESTINATION ${PROJECT_SOURCE_DIR}/lib/}
  COMPONENT library
  )


