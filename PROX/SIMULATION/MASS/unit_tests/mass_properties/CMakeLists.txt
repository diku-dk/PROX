INCLUDE_DIRECTORIES( 
  ${PROJECT_SOURCE_DIR}/PROX/SIMULATION/MASS/MASS/include 
  ${Boost_INCLUDE_DIRS}
)

add_executable(
  unit_mass_properties
  mass_properties.cpp
  )

target_link_libraries(
  unit_mass_properties
  Boost::unit_test_framework
  mass
)

add_test(
  unit_mass_properties
  unit_mass_properties
  )

