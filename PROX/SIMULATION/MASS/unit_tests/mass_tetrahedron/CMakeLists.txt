INCLUDE_DIRECTORIES( 
  ${PROJECT_SOURCE_DIR}/PROX/SIMULATION/MASS/MASS/include 
  ${Boost_INCLUDE_DIRS}
)

add_executable(
  unit_mass_tetrahedron
  mass_tetrahedron.cpp
  )

target_link_libraries(
  unit_mass_tetrahedron
  Boost::unit_test_framework
  mass
)

add_test(
  unit_mass_tetrahedron
  unit_mass_tetrahedron
  )

