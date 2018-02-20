INCLUDE_DIRECTORIES( 
  ${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/GEOMETRY/GEOMETRY/include 
  ${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/TINY/TINY/include 
  ${Boost_INCLUDE_DIRS}
  )

ADD_EXECUTABLE(
  unit_geometry_contacts_obb_cylinder
  geometry_contacts_obb_cylinder.cpp
  )

TARGET_LINK_LIBRARIES(
  unit_geometry_contacts_obb_cylinder
  tiny
  ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
  )

ADD_TEST( 
  unit_geometry_contacts_obb_cylinder
  unit_geometry_contacts_obb_cylinder
  )

