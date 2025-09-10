INCLUDE_DIRECTORIES(
  ${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/GEOMETRY/GEOMETRY/include
${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/EIGENHELPER/EIGENHELPER/include
  ${Boost_INCLUDE_DIRS}
  )

ADD_EXECUTABLE(
  unit_geometry_direction_table
  geometry_direction_table.cpp
  )

TARGET_LINK_LIBRARIES(
  unit_geometry_direction_table
  eigenhelper
  ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
  )

ADD_TEST(
  unit_geometry_direction_table
  unit_geometry_direction_table
  )

