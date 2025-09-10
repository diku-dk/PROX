INCLUDE_DIRECTORIES(
  ${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/GEOMETRY/GEOMETRY/include
${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/EIGENHELPER/EIGENHELPER/include
  ${Boost_INCLUDE_DIRS}
)

ADD_EXECUTABLE(
  unit_geometry_inside_sphere
  geometry_inside_sphere.cpp
  )

TARGET_LINK_LIBRARIES(
  unit_geometry_inside_sphere
  eigenhelper
  geometry
  ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
)

ADD_TEST(
  unit_geometry_inside_sphere
  unit_geometry_inside_sphere
  )

