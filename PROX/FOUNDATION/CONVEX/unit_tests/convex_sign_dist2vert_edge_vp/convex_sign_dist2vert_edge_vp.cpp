#include <tiny.h>
#include <convex.h>

#define BOOST_AUTO_TEST_MAIN
#include <boost/test/auto_unit_test.hpp>
#include <boost/test/unit_test_suite.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/test_tools.hpp>

BOOST_AUTO_TEST_SUITE(convex_signed_distance_to_vertex_edge_vp);

BOOST_AUTO_TEST_CASE(case_by_case_test)
{
  typedef tiny::MathTypes<double>     math_types;
  typedef math_types::vector3_type   vector3_type;
  typedef math_types::real_type      real_type;


  vector3_type a = vector3_type::make(1.0, 0.0, 0.0);
  vector3_type b = vector3_type::make(0.0, 0.0, 0.0);

  // First we use a test point that does not lie on the line

  // Front side of A voronoi plane
  {
    vector3_type p = vector3_type::make( 2.0, 1.0,  1.0);
    real_type sign_p = 0.0; 
    sign_p = convex::signed_distance_to_vertex_edge_voronoi_plane(p, a, b);
    BOOST_CHECK_CLOSE( sign_p, 1.0, 0.01 );
  }
  // Back side of A voronoi plane
  {
    vector3_type p = vector3_type::make( 0.0, 1.0,  1.0);
    real_type sign_p = 0.0; 
    sign_p = convex::signed_distance_to_vertex_edge_voronoi_plane(p, a, b);
    BOOST_CHECK_CLOSE( sign_p, -1.0, 0.01 );
  }
  // In A voronoi plane
  {
    vector3_type p = vector3_type::make( 1.0, 1.0,  1.0);
    real_type sign_p = 0.0; 
    sign_p = convex::signed_distance_to_vertex_edge_voronoi_plane(p, a, b);
    BOOST_CHECK_CLOSE( sign_p, 0.0, 0.01 );
  }

  // Front side of B voronoi plane
  {
    vector3_type p = vector3_type::make( -1.0, 1.0,  1.0);
    real_type sign_p = 0.0; 
    sign_p = convex::signed_distance_to_vertex_edge_voronoi_plane(p, b, a);
    BOOST_CHECK_CLOSE( sign_p, 1.0, 0.01 );
  }
  // Back side of B voronoi plane
  {
    vector3_type p = vector3_type::make( 1.0, 1.0,  1.0);
    real_type sign_p = 0.0; 
    sign_p = convex::signed_distance_to_vertex_edge_voronoi_plane(p, b, a);
    BOOST_CHECK_CLOSE( sign_p, -1.0, 0.01 );
  }
  // In B voronoi plane
  {
    vector3_type p = vector3_type::make( 0.0, 1.0,  1.0);
    real_type sign_p = 0.0; 
    sign_p = convex::signed_distance_to_vertex_edge_voronoi_plane(p, b, a);
    BOOST_CHECK_CLOSE( sign_p, 0.0, 0.01 );
  }

  // Second we use a test point that lies on the line

  // Front side of A voronoi plane
  {
    vector3_type p = vector3_type::make( 2.0, 0.0,  0.0);
    real_type sign_p = 0.0; 
    sign_p = convex::signed_distance_to_vertex_edge_voronoi_plane(p, a, b);
    BOOST_CHECK_CLOSE( sign_p, 1.0, 0.01 );
  }
  // Back side of A voronoi plane
  {
    vector3_type p = vector3_type::make( 0.0, 0.0,  0.0);
    real_type sign_p = 0.0; 
    sign_p = convex::signed_distance_to_vertex_edge_voronoi_plane(p, a, b);
    BOOST_CHECK_CLOSE( sign_p, -1.0, 0.01 );
  }
  // In A voronoi plane
  {
    vector3_type p = vector3_type::make( 1.0, 0.0,  0.0);
    real_type sign_p = 0.0; 
    sign_p = convex::signed_distance_to_vertex_edge_voronoi_plane(p, a, b);
    BOOST_CHECK_CLOSE( sign_p, 0.0, 0.01 );
  }

  // Front side of B voronoi plane
  {
    vector3_type p = vector3_type::make( -1.0, 0.0,  0.0);
    real_type sign_p = 0.0; 
    sign_p = convex::signed_distance_to_vertex_edge_voronoi_plane(p, b, a);
    BOOST_CHECK_CLOSE( sign_p, 1.0, 0.01 );
  }
  // Back side of B voronoi plane
  {
    vector3_type p = vector3_type::make( 1.0, 0.0,  0.0);
    real_type sign_p = 0.0; 
    sign_p = convex::signed_distance_to_vertex_edge_voronoi_plane(p, b, a);
    BOOST_CHECK_CLOSE( sign_p, -1.0, 0.01 );
  }
  // In B voronoi plane
  {
    vector3_type p = vector3_type::make( 0.0, 0.0,  0.0);
    real_type sign_p = 0.0; 
    sign_p = convex::signed_distance_to_vertex_edge_voronoi_plane(p, b, a);
    BOOST_CHECK_CLOSE( sign_p, 0.0, 0.01 );
  }

}

BOOST_AUTO_TEST_SUITE_END();
