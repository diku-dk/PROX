# Try to find atlas (http://math-atlas.sourceforge.net/).
IF(WIN32)

  MESSAGE(ERROR "Unsupported platform?")

ELSE(WIN32)

  FIND_LIBRARY(ATLAS_LIB atlas
    /usr/lib
    /usr/local/lib
    /opt/local/lib
    DOC "What is the path where the file libatlas.la can be found"
    )

  FIND_LIBRARY(LAPACK_LIB lapack
    /usr/lib
    /usr/local/lib
    /opt/local/lib
    DOC "What is the path where the file liblapack.la can be found"
    )

  FIND_LIBRARY(LAPACK_ATLAS_LIB lapack_atlas
    /usr/lib
    /usr/local/lib
    /opt/local/lib
    DOC "What is the path where the file liblapack_atlas.a can be found"
    )

  FIND_LIBRARY(CBLAS_LIB cblas
    /usr/lib
    /usr/local/lib
    /opt/local/lib
    DOC "What is the path where the file libcblas.la can be found"
    )

  FIND_PATH(ATLAS_INCLUDE_DIR atlas/cblas.h
    /usr/include
    /usr/local/include
    /opt/local/include
    DOC "What is the path where the file atlas/cblas.h can be found"
    )

ENDIF(WIN32)

SET(
  ATLAS_LIBRARIES
  ${ATLAS_LIB}
  ${LAPACK_LIB}
  ${CBLAS_LIB}
  ${LAPACK_ATLAS_LIB}
  )

IF(ATLAS_INCLUDE_DIR)
  SET( 
    ATLAS_FOUND  1  
    CACHE STRING
    "Did the system have ATLAS"
    )  
ELSE()
  SET( 
    ATLAS_FOUND  0  
    CACHE STRING 
    "Did the system have ATLAS"
    )
ENDIF()

MARK_AS_ADVANCED( ATLAS_FOUND       )
MARK_AS_ADVANCED( ATLAS_INCLUDE_DIR )
MARK_AS_ADVANCED( ATLAS_LIB         )
MARK_AS_ADVANCED( LAPACK_LIB        )
MARK_AS_ADVANCED( CBLAS_LIB         )
MARK_AS_ADVANCED( LAPACK_ATLAS_LIB  )
MARK_AS_ADVANCED( ATLAS_LIBRARIES   )
