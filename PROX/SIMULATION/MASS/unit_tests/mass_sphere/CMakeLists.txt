INCLUDE_DIRECTORIES( 
  ${PROJECT_SOURCE_DIR}/PROX/SIMULATION/MASS/MASS/include 
  ${Boost_INCLUDE_DIRS}
)

add_executable(
  unit_mass_sphere
  mass_sphere.cpp
  )

target_link_libraries(
  unit_mass_sphere 
  Boost::unit_test_framework
  mass
)

add_test(
  unit_mass_sphere
  unit_mass_sphere
  )

