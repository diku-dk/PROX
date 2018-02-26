INCLUDE_DIRECTORIES( 
  ${PROJECT_SOURCE_DIR}/FOUNDATION/SPARSE/SPARSE/include 
  ${PROJECT_SOURCE_DIR}/FOUNDATION/UTIL/UTIL/include 
  ${Boost_INCLUDE_DIRS}
  )

ADD_EXECUTABLE(
  unit_sparse_fake_transpose
  sparse_fake_transpose.cpp
  )

TARGET_LINK_LIBRARIES(
  unit_sparse_fake_transpose
  util
  ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
  )

ADD_TEST(
  unit_sparse_fake_transpose
  unit_sparse_fake_transpose
  )