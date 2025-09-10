INCLUDE_DIRECTORIES(
  ${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/CONVEX/CONVEX/include
${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/EIGENHELPER/EIGENHELPER/include
  ${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/GEOMETRY/GEOMETRY/include
  ${Boost_INCLUDE_DIRS}
)

ADD_EXECUTABLE(
  unit_convex_conservative_advancement
  convex_conservative_advancement.cpp
  )

TARGET_LINK_LIBRARIES(
  unit_convex_conservative_advancement
  eigenhelper
  geometry
  convex
  ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
  )

ADD_TEST(
  unit_convex_conservative_advancement
  unit_convex_conservative_advancement
  )

