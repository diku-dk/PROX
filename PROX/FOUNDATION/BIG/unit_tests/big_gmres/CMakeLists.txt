INCLUDE_DIRECTORIES( 
  ${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/BIG/BIG/include
  ${Boost_INCLUDE_DIRS} 
  )

ADD_EXECUTABLE(
  unit_big_gmres 
  unit_big_gmres.cpp
  )

TARGET_LINK_LIBRARIES(
  unit_big_gmres
  big
  ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
)

ADD_TEST( 
  unit_big_gmres 
  unit_big_gmres 
  )

