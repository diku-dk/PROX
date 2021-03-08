add_executable(
  unit_mass_transform
  mass_transform.cpp
  )

target_link_libraries(
  unit_mass_transform
  Boost::unit_test_framework
  tiny
  mass
)

add_test(
  unit_mass_transform
  unit_mass_transform
  )

