INCLUDE_DIRECTORIES( 
  ${PROJECT_SOURCE_DIR}/FOUNDATION/UTIL/UTIL/include 
  ${PROJECT_SOURCE_DIR}/FOUNDATION/TINY/TINY/include 
  ${PROJECT_SOURCE_DIR}/FOUNDATION/GEOMETRY/GEOMETRY/include 
  ${PROJECT_SOURCE_DIR}/FOUNDATION/MESH_ARRAY/MESH_ARRAY/include 
  ${PROJECT_SOURCE_DIR}/SIMULATION/KDOP/KDOP/include 
  ${Boost_INCLUDE_DIRS}
  )

ADD_EXECUTABLE(
  unit_kdop_traversal
  kdop_traversal.cpp
  )

TARGET_LINK_LIBRARIES(
  unit_kdop_traversal
  tiny
  geometry
  mesh_array
  tetgen
  util
  ${Boost_SYSTEM_LIBRARY}
  ${Boost_FILESYSTEM_LIBRARY}
  ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
  )

ADD_TEST( 
  unit_kdop_traversal
  unit_kdop_traversal
  )
