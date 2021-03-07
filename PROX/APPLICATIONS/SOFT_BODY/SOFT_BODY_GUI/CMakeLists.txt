set( CMAKE_INCLUDE_CURRENT_DIR ON                                     )
set( CMAKE_AUTOMOC ON                                                 )

set( SOFT_BODY_GUI_CODE_PATH     "${PROJECT_SOURCE_DIR}/PROX/APPLICATIONS/SOFT_BODY/SOFT_BODY_GUI"          )
set( SOFT_BODY_GUI_NAME          "soft_body_gui"                                                       )
set( SOFT_BODY_GUI_RUN_PATH      "${PROJECT_SOURCE_DIR}/bin/soft_body_gui"                             )
set( SOFT_BODY_GUI_SHADERS_PATH  "${PROJECT_SOURCE_DIR}/PROX/APPLICATIONS/SOFT_BODY/SOFT_BODY_GUI/shaders"  )
set( SOFT_BODY_GUI_TEXTURES_PATH "${PROJECT_SOURCE_DIR}/bin/resources/textures"                            )
set( SOFT_BODY_GUI_OBJS_PATH     "${PROJECT_SOURCE_DIR}/bin/resources/objs"                                )
set( SOFT_BODY_GUI_CFGS_PATH     "${PROJECT_SOURCE_DIR}/bin/soft_body_gui"                             )

file(GLOB 
  SOFT_BODY_GUI_SOURCES
  ${SOFT_BODY_GUI_CODE_PATH}/src/*.cpp 
  )

file(GLOB
  SOFT_BODY_GUI_HEADERS
  ${SOFT_BODY_GUI_CODE_PATH}/include/*.h 
  )

file(GLOB
  SOFT_BODY_GUI_SHADERS
  ${SOFT_BODY_GUI_SHADERS_PATH}/*.glsl
  )

file(GLOB
  SOFT_BODY_GUI_CFGS
  ${SOFT_BODY_GUI_CFGS_PATH}/*.cfg
  )

add_executable(${SOFT_BODY_GUI_NAME} 
  ${SOFT_BODY_GUI_SOURCES} 
  ${SOFT_BODY_GUI_HEADERS}
  ${SOFT_BODY_GUI_SHADERS}
  ${SOFT_BODY_GUI_CFGS}
  )
  
source_group(Shaders FILES ${SOFT_BODY_GUI_SHADERS} )

source_group(Config  FILES ${SOFT_BODY_GUI_CFGS}    )

target_include_directories(${SOFT_BODY_GUI_NAME}
  PRIVATE
  ${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/GL3/GL3/include    
  ${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/UTIL/UTIL/include
  ${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/TINY/TINY/include 
  ${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/SPARSE/SPARSE/include 
  ${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/GEOMETRY/GEOMETRY/include 
  ${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/MESH_ARRAY/MESH_ARRAY/include
  ${PROJECT_SOURCE_DIR}/PROX/SIMULATION/KDOP/KDOP/include
  ${PROJECT_SOURCE_DIR}/PROX/SIMULATION/HYPER/HYPER/include    
  ${PROJECT_SOURCE_DIR}/PROX/APPLICATIONS/SOFT_BODY/SOFT_BODY_SHARED/include 
  ${SOFT_BODY_GUI_CODE_PATH}/include
  )

target_link_libraries(${SOFT_BODY_GUI_NAME}
  PRIVATE
  soft_body_shared
  gl3
  kdop
  hyper
  tiny
  sparse
  geometry
  mesh_array
  util
  tetgen
  qhull
  Qt5::Widgets
  Qt5::Gui
  Qt5::Core
  Qt5::OpenGL
  OpenGL::GL
  GLEW::glew
  glm::glm
  Boost::system
  Boost::filesystem
  )

target_compile_definitions(${SOFT_BODY_GUI_NAME}
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
  TARGETS ${SOFT_BODY_GUI_NAME}
  RUNTIME DESTINATION ${PROJECT_SOURCE_DIR}/bin/}
  COMPONENT binary
  )
