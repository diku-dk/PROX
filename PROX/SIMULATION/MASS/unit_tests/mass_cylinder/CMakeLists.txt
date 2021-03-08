INCLUDE_DIRECTORIES( 
  ${PROJECT_SOURCE_DIR}/PROX/SIMULATION/MASS/MASS/include 
  ${Boost_INCLUDE_DIRS}
)

add_executable(
  unit_mass_cylinder
  mass_cylinder.cpp
  )

target_link_libraries(
  unit_mass_cylinder
  Boost::unit_test_framework
  mass
)

add_test(
  unit_mass_cylinder
  unit_mass_cylinder
  )

