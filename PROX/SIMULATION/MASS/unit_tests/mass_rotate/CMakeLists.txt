INCLUDE_DIRECTORIES( 
  ${PROJECT_SOURCE_DIR}/PROX/FOUNDATION/TINY/TINY/include 
  ${PROJECT_SOURCE_DIR}/PROX/SIMULATION/MASS/MASS/include 
  ${Boost_INCLUDE_DIRS}
)

add_executable(
  unit_mass_rotate
  mass_rotate.cpp
  )

target_link_libraries(
  unit_mass_rotate
  Boost::unit_test_framework
  tiny
  mass
)

add_test(
  unit_mass_rotate
  unit_mass_rotate
  )

