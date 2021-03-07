add_library(tinyxml
      tinyxml.cpp
      tinyxmlerror.cpp
      tinyxmlparser.cpp
	    tinyxml.h
	    tinystr.h
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
  RUNTIME DESTINATION ${PROJECT_SOURCE_DIR}/lib/}
  COMPONENT library
 )

