INCLUDE_DIRECTORIES( 
  ${PROJECT_SOURCE_DIR}/FOUNDATION/BIG/BIG/include
  ${Boost_INCLUDE_DIRS} 
  )

ADD_EXECUTABLE(
  unit_big_gram_schmidt 
  unit_big_gram_schmidt.cpp
  )

TARGET_LINK_LIBRARIES(
  unit_big_gram_schmidt
  big
  ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
)

ADD_TEST( 
  unit_big_gram_schmidt 
  unit_big_gram_schmidt 
  )


