INCLUDE_DIRECTORIES( 
  ${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/UTIL/UTIL/include 
  ${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/TINY/TINY/include 
  ${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/GEOMETRY/GEOMETRY/include 
  ${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/MESH_ARRAY/MESH_ARRAY/include 
  ${PROJECT_SOURCE_DIR}/PROX/SIMULATION/KDOP/KDOP/include 
  ${Boost_INCLUDE_DIRS}
  )

add_executable(
  unit_kdop_raycast
  kdop_raycast.cpp
  )

target_link_libraries(
  unit_kdop_raycast
  tiny
  geometry
  mesh_array
  tetgen
  util
  ${Boost_SYSTEM_LIBRARY}
  ${Boost_FILESYSTEM_LIBRARY}
  Boost::unit_test_framework
  )

add_test( 
  unit_kdop_raycast
  unit_kdop_raycast
  )
