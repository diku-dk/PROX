add_library(tinyxml
      tinyxml.cpp
      tinyxmlerror.cpp
      tinyxmlparser.cpp
	    tinyxml.h
	    tinystr.h
      )

set_target_properties(tinyxml PROPERTIES PUBLIC_HEADER "tinyxml.h;tinystr.h")

target_include_directories(tinyxml
  PUBLIC
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/PROX/3RDPARTY/TinyXML>
  $<INSTALL_INTERFACE:include>
  )

#-------------------------------------------------------
#
# Tell that TinyXML in 3RDPARTY should use
# std::strings over char *
#
# TinyXML can work both with stl strings and char*. The
# functionality is controlled by defining the pre-processor
# directive TIXML_USE_STL or not.
#
# Note that every target that links against TinyXML
# also need to have the TIXML_USE_STL flag defined!
#

target_compile_definitions(tinyxml
  PUBLIC
  TIXML_USE_STL
)

install(
  TARGETS tinyxml
  ARCHIVE DESTINATION lib
  PUBLIC_HEADER DESTINATION include
 )

