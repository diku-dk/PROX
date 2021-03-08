INCLUDE_DIRECTORIES( 
  ${PROJECT_SOURCE_DIR}/PROX/SIMULATION/MASS/MASS/include 
  ${Boost_INCLUDE_DIRS}
)

add_executable(
  unit_mass_translate
  mass_translate.cpp
  )

target_link_libraries(
  unit_mass_translate
  Boost::unit_test_framework
  mass
)

add_test(
  unit_mass_translate
  unit_mass_translate
  )

