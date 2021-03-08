file(GLOB 
  RIGID_BODY_GUI_SHARED_SOURCES
  src/*.cpp 
  )

file(GLOB
  RIGID_BODY_GUI_SHARED_HEADERS
  include/*.h 
  )

file(GLOB
  RIGID_BODY_GUI_SHARED_SHADERS
  shaders/*.glsl
  shaders/ssao/*.glsl
  shaders/show_texture/*.glsl
  )

add_library(rigid_body_gui_shared 
  ${RIGID_BODY_GUI_SHARED_SOURCES} 
  ${RIGID_BODY_GUI_SHARED_HEADERS}
  ${RIGID_BODY_GUI_SHARED_SHADERS}
  )

source_group(Shaders FILES ${RIGID_BODY_GUI_SHARED_SHADERS} )

target_include_directories(rigid_body_gui_shared
  PUBLIC
  ${PROJECT_SOURCE_DIR}/PROX/APPLICATIONS/RIGID_BODY/RIGID_BODY_GUI_SHARED/include
)

target_link_libraries(rigid_body_gui_shared
  PUBLIC
  gl3
  content
  tiny
  util
  mesh_array
)


target_compile_definitions(rigid_body_gui_shared
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
  TARGETS rigid_body_gui_shared
  RUNTIME DESTINATION  ${RIGID_BODY_GUI_RUN_PATH}
  COMPONENT binary
 )
