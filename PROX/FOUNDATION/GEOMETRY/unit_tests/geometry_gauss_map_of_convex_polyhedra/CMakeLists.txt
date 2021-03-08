add_executable(
  unit_geometry_gauss_map_of_convex_polyhedra
  geometry_gauss_map_of_convex_polyhedra.cpp
  )

target_link_libraries(
  unit_geometry_gauss_map_of_convex_polyhedra
  geometry
  Boost::unit_test_framework
  )

add_test( 
  unit_geometry_gauss_map_of_convex_polyhedra
  unit_geometry_gauss_map_of_convex_polyhedra
  )
