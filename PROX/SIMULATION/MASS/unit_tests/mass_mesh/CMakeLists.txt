add_executable(
  unit_mass_mesh
  mass_mesh.cpp
  )

target_link_libraries(
  unit_mass_mesh
  mass
  Boost::unit_test_framework
  )

add_test(
  unit_mass_mesh
  unit_mass_mesh
  )
