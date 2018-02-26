INCLUDE_DIRECTORIES( 
  ${PROJECT_SOURCE_DIR}/FOUNDATION/BIG/BIG/include
  ${Boost_INCLUDE_DIRS} 
  )

ADD_EXECUTABLE(
  unit_big_prod_row 
  unit_big_prod_row.cpp
  )

TARGET_LINK_LIBRARIES(
  unit_big_prod_row
  big
  ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
)

ADD_TEST( 
  unit_big_prod_row 
  unit_big_prod_row 
  )

