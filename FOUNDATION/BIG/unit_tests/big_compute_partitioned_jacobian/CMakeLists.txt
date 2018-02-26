INCLUDE_DIRECTORIES( 
  ${PROJECT_SOURCE_DIR}/FOUNDATION/BIG/BIG/include
  ${Boost_INCLUDE_DIRS} 
  )

ADD_EXECUTABLE(
  unit_big_compute_partitioned_jacobian 
  unit_big_compute_partitioned_jacobian.cpp
  )

TARGET_LINK_LIBRARIES(
  unit_big_compute_partitioned_jacobian
  big
  ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
)

ADD_TEST( 
  unit_big_compute_partitioned_jacobian 
  unit_big_compute_partitioned_jacobian 
  )

