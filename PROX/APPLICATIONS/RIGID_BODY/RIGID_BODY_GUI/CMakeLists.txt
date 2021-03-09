set( CMAKE_INCLUDE_CURRENT_DIR ON )
set( CMAKE_AUTOMOC ON )

file(GLOB 
  RIGID_BODY_GUI_SOURCES
  ${PROJECT_SOURCE_DIR}/PROX/APPLICATIONS/RIGID_BODY/RIGID_BODY_GUI/src/*.cpp 
  )

file(GLOB
  RIGID_BODY_GUI_HEADERS
  ${PROJECT_SOURCE_DIR}/PROX/APPLICATIONS/RIGID_BODY/RIGID_BODY_GUI/include/*.h 
  )

file(GLOB
  RIGID_BODY_GUI_CFGS
  ${PROJECT_SOURCE_DIR}/bin/rigid_body_gui/*.cfg
  )

add_executable(rigid_body_gui
  ${RIGID_BODY_GUI_SOURCES} 
  ${RIGID_BODY_GUI_HEADERS}
  ${RIGID_BODY_GUI_CFGS}
  )

source_group(Config FILES ${RIGID_BODY_GUI_CFGS})

target_include_directories(rigid_body_gui
  PRIVATE
  ${PROJECT_SOURCE_DIR}/PROX/APPLICATIONS/RIGID_BODY/RIGID_BODY_GUI/include
)

target_link_libraries(rigid_body_gui
  PRIVATE
  rigid_body_gui_shared
  gl3
  tiny
  simulators
  procedural
  util
  Qt5::Widgets
  Qt5::Gui
  Qt5::Core
  Qt5::OpenGL
  Boost::system
  Boost::filesystem
)

target_compile_definitions(rigid_body_gui
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
  TARGETS rigid_body_gui
  RUNTIME DESTINATION bin
 )
