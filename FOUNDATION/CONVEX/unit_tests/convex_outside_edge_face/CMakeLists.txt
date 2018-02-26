INCLUDE_DIRECTORIES( 
  ${PROJECT_SOURCE_DIR}/FOUNDATION/CONVEX/CONVEX/include 
  ${PROJECT_SOURCE_DIR}/FOUNDATION/CONVEX/CONVEX/private 
  ${PROJECT_SOURCE_DIR}/FOUNDATION/TINY/TINY/include 
  ${PROJECT_SOURCE_DIR}/FOUNDATION/GEOMETRY/GEOMETRY/include 
  ${Boost_INCLUDE_DIRS}
)

ADD_EXECUTABLE(
  unit_convex_outside_edge_face
  convex_outside_edge_face.cpp
  )

TARGET_LINK_LIBRARIES(
  unit_convex_outside_edge_face 
  tiny
  geometry
  convex
  ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
)

ADD_TEST(
  unit_convex_outside_edge_face
  unit_convex_outside_edge_face
  )

