INCLUDE_DIRECTORIES( 
  ${PROJECT_SOURCE_DIR}/PROX/SIMULATION/MASS/MASS/include 
  ${Boost_INCLUDE_DIRS}
)

add_executable(
  unit_mass_box
  mass_box.cpp
  )

target_link_libraries(
  unit_mass_box 
  Boost::unit_test_framework
  mass
)

add_test(
  unit_mass_box
  unit_mass_box
  )

