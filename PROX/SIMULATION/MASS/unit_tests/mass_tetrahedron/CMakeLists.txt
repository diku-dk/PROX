add_executable(
  unit_mass_tetrahedron
  mass_tetrahedron.cpp
  )

target_link_libraries(
  unit_mass_tetrahedron
  mass
  Boost::unit_test_framework
  )

add_test(
  unit_mass_tetrahedron
  unit_mass_tetrahedron
  )
