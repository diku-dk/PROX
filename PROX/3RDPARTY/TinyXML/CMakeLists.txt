add_library(tinyxml
            tinyxml.cpp  
            tinyxmlerror.cpp  
            tinyxmlparser.cpp
	    tinyxml.h
	    tinystr.h
            )

install(
  TARGETS tinyxml
  RUNTIME DESTINATION ${PROJECT_SOURCE_DIR}/lib/}
  COMPONENT library
 )

