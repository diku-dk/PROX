set( CMAKE_INCLUDE_CURRENT_DIR ON                                     )
set( CMAKE_AUTOMOC ON                                                 )

file(GLOB 
  SOFT_BODY_GUI_SOURCES
  ${PROJECT_SOURCE_DIR}/PROX/APPLICATIONS/SOFT_BODY/SOFT_BODY_GUI/src/*.cpp
  )

file(GLOB
  SOFT_BODY_GUI_HEADERS
  ${PROJECT_SOURCE_DIR}/PROX/APPLICATIONS/SOFT_BODY/SOFT_BODY_GUI/include/*.h
  )

file(GLOB
  SOFT_BODY_GUI_SHADERS
  ${PROJECT_SOURCE_DIR}/PROX/APPLICATIONS/SOFT_BODY/SOFT_BODY_GUI/shaders/*.glsl
  )

file(GLOB
  SOFT_BODY_GUI_CFGS
  ${PROJECT_SOURCE_DIR}/bin/soft_body_gui/*.cfg
  )

add_executable(soft_body_gui
  ${SOFT_BODY_GUI_SOURCES} 
  ${SOFT_BODY_GUI_HEADERS}
  ${SOFT_BODY_GUI_SHADERS}
  ${SOFT_BODY_GUI_CFGS}
  )
  
source_group(Shaders FILES ${SOFT_BODY_GUI_SHADERS} )

source_group(Config  FILES ${SOFT_BODY_GUI_CFGS}    )

target_include_directories(soft_body_gui
  PRIVATE
  ${PROJECT_SOURCE_DIR}/PROX/APPLICATIONS/SOFT_BODY/SOFT_BODY_GUI/include
  )

target_link_libraries(soft_body_gui
  PRIVATE
  soft_body_shared
  gl3
  tiny
  geometry
  mesh_array
  hyper
  util
  Qt5::Widgets
  Qt5::Gui
  Qt5::Core
  Qt5::OpenGL
  Boost::system
  Boost::filesystem
  )

target_compile_definitions(soft_body_gui
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
  TARGETS soft_body_gui
  RUNTIME DESTINATION ${PROJECT_SOURCE_DIR}/bin/}
  COMPONENT binary
  )
