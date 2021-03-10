file(GLOB 
  SOFT_BODY_SHARED_SOURCES
  src/*.cpp 
  )

file(GLOB
  SOFT_BODY_SHARED_HEADERS
  include/*.h 
  )

add_library(soft_body_shared
  ${SOFT_BODY_SHARED_SOURCES} 
  ${SOFT_BODY_SHARED_HEADERS}
  )

set_target_properties(soft_body_shared PROPERTIES PUBLIC_HEADER "${SOFT_BODY_SHARED_HEADERS}")

target_include_directories(soft_body_shared
  PUBLIC
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/PROX/APPLICATIONS/SOFT_BODY/SOFT_BODY_SHARED/include>
  $<INSTALL_INTERFACE:include>
  )

target_link_libraries(soft_body_shared
  PUBLIC
  hyper
  tiny
  util
  mesh_array
  )

target_compile_definitions(soft_body_shared
  PUBLIC
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
  TARGETS soft_body_shared
  ARCHIVE DESTINATION lib
  PUBLIC_HEADER DESTINATION include
  )
