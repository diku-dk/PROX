#include <mesh_array.h>
#include <tiny.h>

#define BOOST_AUTO_TEST_MAIN
#include <boost/test/auto_unit_test.hpp>
#include <boost/test/unit_test_suite.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/test_tools.hpp>

#include <cmath>   // needed for std::sqrt and std::fabs

BOOST_AUTO_TEST_SUITE(mesh_array);

BOOST_AUTO_TEST_CASE(mesh_array_compute_distance_map_sphere)
{
  using std::sqrt;
  using std::fabs;

  typedef tiny::MathTypes<float> MT;
  typedef MT::real_type          T;
  
  mesh_array::T3Mesh surf;
  mesh_array::VertexAttribute<T,mesh_array::T3Mesh> sX;
  mesh_array::VertexAttribute<T,mesh_array::T3Mesh> sY;
  mesh_array::VertexAttribute<T,mesh_array::T3Mesh> sZ;

  size_t const slices   = 12u;
  size_t const segments = 24u;
  T      const radius   = 12.0f;
  
  mesh_array::make_sphere<MT>(radius, slices, segments, surf, sX, sY, sZ);
  
  mesh_array::T4Mesh mesh;
  mesh_array::VertexAttribute<T,mesh_array::T4Mesh> X;
  mesh_array::VertexAttribute<T,mesh_array::T4Mesh> Y;
  mesh_array::VertexAttribute<T,mesh_array::T4Mesh> Z;

  mesh_array::tetgen(surf, sX, sY, sZ, mesh, X, Y, Z, mesh_array::tetgen_quality_settings() );
  
  mesh_array::VertexAttribute<T,mesh_array::T4Mesh> phi;
  
  mesh_array::compute_distance_map( mesh, X, Y, Z, phi );
  
  for (size_t i = 0u; i < phi.size(); ++i)
  {
    T const real_distance = radius - std::sqrt(X(i)*X(i) + Y(i)*Y(i) + Z(i)*Z(i));

    BOOST_CHECK_CLOSE(  phi(i), real_distance, 1.0 );
  }
}

BOOST_AUTO_TEST_SUITE_END();
