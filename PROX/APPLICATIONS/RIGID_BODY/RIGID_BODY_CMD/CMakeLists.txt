file(GLOB
  RIGID_BODY_CMD_SOURCES
  src/*.cpp 
  )

file(GLOB
  RIGID_BODY_CMD_HEADERS
  include/*.h 
  )

file(GLOB
  RIGID_BODY_CMD_CFGS
  ${PROJECT_SOURCE_DIR}/bin/rigid_body_cmd/*.cfg
  )

add_executable(rigid_body_cmd 
  ${RIGID_BODY_CMD_SOURCES} 
  ${RIGID_BODY_CMD_HEADERS}
  ${RIGID_BODY_CMD_CFGS}
  )

source_group(Config FILES ${RIGID_BODY_CMD_CFGS})

target_include_directories(rigid_body_cmd
  PRIVATE
  ${PROJECT_SOURCE_DIR}/PROX/APPLICATIONS/RIGID_BODY/RIGID_BODY_CMD/include 
  )

target_link_libraries(rigid_body_cmd
  PRIVATE
  tiny
  mass
  geometry
  mesh_array
  simulators
  procedural
  util
  Boost::boost
  )

target_compile_definitions(rigid_body_cmd
  PRIVATE
    $<$<CXX_COMPILER_ID:MSVC>:_SCL_SECURE_NO_WARNINGS>
    $<$<CXX_COMPILER_ID:MSVC>:_SCL_SECURE_NO_DEPRECATE>
    $<$<CXX_COMPILER_ID:MSVC>:_CRT_SECURE_NO_WARNINGS>
    $<$<CXX_COMPILER_ID:MSVC>:_CRT_SECURE_NO_DEPRECATE>
    $<$<CXX_COMPILER_ID:MSVC>:_SILENCE_STDEXT_HASH_DEPRECATION_WARNINGS>
    $<$<CXX_COMPILER_ID:GNU>:_USE_MATH_DEFINES>
    $<$<CXX_COMPILER_ID:GNU>:GL_SILENCE_DEPRECATION>
    $<$<CXX_COMPILER_ID:AppleClang>:_USE_MATH_DEFINES>
    $<$<CXX_COMPILER_ID:AppleClang>:GL_SILENCE_DEPRECATION>
)

install(
  TARGETS rigid_body_cmd
  RUNTIME DESTINATION bin
  )
