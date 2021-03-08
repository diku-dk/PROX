INCLUDE_DIRECTORIES( 
  ${PROJECT_SOURCE_DIR}/PROX/SIMULATION/MASS/MASS/include 
  ${Boost_INCLUDE_DIRS}
)

add_executable(
  unit_mass_mesh
  mass_mesh.cpp
  )

target_link_libraries(
  unit_mass_mesh
  Boost::unit_test_framework
  mass
)

add_test(
  unit_mass_mesh
  unit_mass_mesh
  )

