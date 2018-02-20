INCLUDE_DIRECTORIES( 
  ${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/BIG/BIG/include
  ${Boost_INCLUDE_DIRS} 
  )

ADD_EXECUTABLE(
  unit_big_prod_sub_rhs 
  unit_big_prod_sub_rhs.cpp
  )

TARGET_LINK_LIBRARIES(
  unit_big_prod_sub_rhs
  big
  ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
)

ADD_TEST( 
  unit_big_prod_sub_rhs 
  unit_big_prod_sub_rhs 
  )
