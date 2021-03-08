INCLUDE_DIRECTORIES( 
  ${PROJECT_SOURCE_DIR}/PROX/SIMULATION/MASS/MASS/include 
  ${Boost_INCLUDE_DIRS}
)

add_executable(
  unit_mass_capsule
  mass_capsule.cpp
  )

target_link_libraries(
  unit_mass_capsule 
  Boost::unit_test_framework
  mass
)

add_test(
  unit_mass_capsule
  unit_mass_capsule
  )

