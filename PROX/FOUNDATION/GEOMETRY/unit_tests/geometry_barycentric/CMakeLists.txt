INCLUDE_DIRECTORIES(
  ${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/GEOMETRY/GEOMETRY/include
  ${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/EIGENHELPER/EIGENHELPER/include
  ${Boost_INCLUDE_DIRS}
  )

ADD_EXECUTABLE(
  unit_geometry_barycentric
  geometry_barycentric.cpp
  )

TARGET_LINK_LIBRARIES(
  unit_geometry_barycentric
  eigenhelper
  ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
  )

ADD_TEST(
  unit_geometry_barycentric
  unit_geometry_barycentric
  )

