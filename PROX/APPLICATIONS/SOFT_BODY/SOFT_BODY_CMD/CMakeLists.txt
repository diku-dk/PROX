file(GLOB 
  SOFT_BODY_CMD_SOURCES
  src/*.cpp 
  )

file(GLOB
  SOFT_BODY_CMD_HEADERS
  include/*.h 
  )
  
file(GLOB
  SOFT_BODY_CMD_CFGS
  ${PROJECT_SOURCE_DIR}/bin/soft_body_cmd/*.cfg
  )

add_executable(soft_body_cmd
  ${SOFT_BODY_CMD_SOURCES} 
  ${SOFT_BODY_CMD_HEADERS}
  ${SOFT_BODY_CMD_CFGS}
  )

source_group(Config    FILES ${SOFT_BODY_CMD_CFGS})

target_include_directories(soft_body_cmd
  PRIVATE
  ${PROJECT_SOURCE_DIR}/PROX/APPLICATIONS/SOFT_BODY/SOFT_BODY_CMD/include
  )

target_link_libraries(soft_body_cmd
  PRIVATE
  soft_body_shared
  hyper
  util
  )

target_compile_definitions(soft_body_cmd
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
  TARGETS soft_body_cmd
  RUNTIME DESTINATION bin
  )

set_target_properties(soft_body_cmd
  PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY_DEBUG   "${PROJECT_SOURCE_DIR}/bin/soft_body_cmd/"
  RUNTIME_OUTPUT_DIRECTORY_RELEASE "${PROJECT_SOURCE_DIR}/bin/soft_body_cmd/"
  )
