INCLUDE_DIRECTORIES( 
  ${PROJECT_SOURCE_DIR}/PROX/SIMULATION/MASS/MASS/include 
  ${Boost_INCLUDE_DIRS}
)

add_executable(
  unit_mass_cone
  mass_cone.cpp
  )

target_link_libraries(
  unit_mass_cone
  Boost::unit_test_framework
  mass
)

add_test(
  unit_mass_cone
  unit_mass_cone
  )

