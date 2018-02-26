INCLUDE_DIRECTORIES( 
  ${PROJECT_SOURCE_DIR}/FOUNDATION/GEOMETRY/GEOMETRY/include 
  ${PROJECT_SOURCE_DIR}/FOUNDATION/TINY/TINY/include 
  ${Boost_INCLUDE_DIRS}
  )

ADD_EXECUTABLE(
  unit_geometry_contacts_obb_obb
  geometry_contacts_obb_obb.cpp
  )

TARGET_LINK_LIBRARIES(
  unit_geometry_contacts_obb_obb
  tiny
  ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
  )

ADD_TEST( 
  unit_geometry_contacts_obb_obb
  unit_geometry_contacts_obb_obb
  )

