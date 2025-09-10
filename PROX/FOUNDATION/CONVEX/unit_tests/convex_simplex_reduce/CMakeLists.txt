INCLUDE_DIRECTORIES(
  ${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/CONVEX/CONVEX/include
  ${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/CONVEX/CONVEX/private
  ${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/EIGENHELPER/EIGENHELPER/include
  ${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/GEOMETRY/GEOMETRY/include
  ${Boost_INCLUDE_DIRS}
)

ADD_EXECUTABLE(
  unit_convex_simplex_reduce
  convex_simplex_reduce.cpp
  )

TARGET_LINK_LIBRARIES(
  unit_convex_simplex_reduce
  eigenhelper
  geometry
  convex
  ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
)

ADD_TEST(
  unit_convex_simplex_reduce
  unit_convex_simplex_reduce
  )

