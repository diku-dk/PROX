INCLUDE_DIRECTORIES(
  ${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/GEOMETRY/GEOMETRY/include
${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/EIGENHELPER/EIGENHELPER/include
  ${Boost_INCLUDE_DIRS}
  )

ADD_EXECUTABLE(
  unit_geometry_overlap_triangle_triangle
  geometry_overlap_triangle_triangle.cpp
  )

TARGET_LINK_LIBRARIES(
  unit_geometry_overlap_triangle_triangle
  eigenhelper
  ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
  )

ADD_TEST(
  unit_geometry_overlap_triangle_triangle
  unit_geometry_overlap_triangle_triangle
  )

