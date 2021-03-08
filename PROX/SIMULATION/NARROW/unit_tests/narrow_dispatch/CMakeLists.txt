INCLUDE_DIRECTORIES( 
  ${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/UTIL/UTIL/include 
  ${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/TINY/TINY/include 
  ${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/GEOMETRY/GEOMETRY/include 
  ${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/MESH_ARRAY/MESH_ARRAY/include
  ${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/CONVEX/CONVEX/include 
  ${PROJECT_SOURCE_DIR}/PROX/SIMULATION/NARROW/NARROW/include
  ${PROJECT_SOURCE_DIR}/PROX/SIMULATION/KDOP/KDOP/include
  ${Boost_INCLUDE_DIRS} 
  )

IF(HAS_DIKUCL)
  INCLUDE_DIRECTORIES(
      ${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/DIKUCL/DIKUCL/include 
    )
ENDIF(HAS_DIKUCL)

add_executable(
  unit_narrow_dispatch
  narrow_dispatch.cpp
  )

target_link_libraries(
  unit_narrow_dispatch
  narrow
  kdop
  tiny
  geometry
  mesh_array
  convex
  ${Boost_SYSTEM_LIBRARY}
  ${Boost_FILESYSTEM_LIBRARY}
  Boost::unit_test_framework
  qhull
  tetgen
  util
  )

IF(HAS_DIKUCL)
  target_link_libraries(unit_narrow_dispatch dikucl)
ENDIF(HAS_DIKUCL)


add_test(
  unit_narrow_dispatch
  unit_narrow_dispatch
  )

