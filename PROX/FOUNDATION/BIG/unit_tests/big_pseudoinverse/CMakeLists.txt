INCLUDE_DIRECTORIES( 
  ${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/BIG/BIG/include
  ${Boost_INCLUDE_DIRS} 
  )

ADD_EXECUTABLE(
  unit_big_pseudoinverse 
  unit_big_pseudoinverse.cpp
  )

TARGET_LINK_LIBRARIES(
  unit_big_pseudoinverse
  big
  ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
)

ADD_TEST( 
  unit_big_pseudoinverse 
  unit_big_pseudoinverse 
  )

