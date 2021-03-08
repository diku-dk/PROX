add_executable(
  unit_mass_ellipsoid
  mass_ellipsoid.cpp
  )

target_link_libraries(
  unit_mass_ellipsoid
  mass
  Boost::unit_test_framework
  )

add_test(
  unit_mass_ellipsoid
  unit_mass_ellipsoid
  )
