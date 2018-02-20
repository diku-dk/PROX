INCLUDE_DIRECTORIES( 
  ${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/CONVEX/CONVEX/include 
  ${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/TINY/TINY/include 
  ${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/GEOMETRY/GEOMETRY/include 
  ${Boost_INCLUDE_DIRS}
)

ADD_EXECUTABLE(
  unit_convex_motion_interpolation
  convex_motion_interpolation.cpp
  )

TARGET_LINK_LIBRARIES(
  unit_convex_motion_interpolation 
  tiny
  geometry
  convex
  ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
  )

ADD_TEST(
  unit_convex_motion_interpolation
  unit_convex_motion_interpolation
  )

