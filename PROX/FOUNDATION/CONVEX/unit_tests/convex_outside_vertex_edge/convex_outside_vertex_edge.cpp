#include <tiny.h>
#include <convex.h>

#define BOOST_AUTO_TEST_MAIN
#include <boost/test/auto_unit_test.hpp>
#include <boost/test/unit_test_suite.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/test_tools.hpp>

BOOST_AUTO_TEST_SUITE(convex_outside_vertex_edge);

BOOST_AUTO_TEST_CASE(case_by_case_test)
{
  typedef tiny::MathTypes<double>     math_types;
  typedef math_types::vector3_type   V;

  V const a = V::make(1.0, 0.0, 0.0);
  V const b = V::make(0.0, 0.0, 0.0);

  // First we use a test point that does not lie on the line

  // Front side of A voronoi plane
  {
    V p = V::make( 2.0, 1.0,  1.0);
    bool outside = convex::outside_vertex_edge_voronoi_plane(p, a, b);
    BOOST_CHECK( outside );
  }
  // Back side of A voronoi plane
  {
    V p = V::make( 0.0, 1.0,  1.0);
    bool outside = convex::outside_vertex_edge_voronoi_plane(p, a, b);
    BOOST_CHECK( !outside );
  }
  // In A voronoi plane
  {
    V p = V::make( 1.0, 1.0,  1.0);
    bool outside = convex::outside_vertex_edge_voronoi_plane(p, a, b);
    BOOST_CHECK( outside );
  }

  // Front side of B voronoi plane
  {
    V p = V::make( -1.0, 1.0,  1.0);
    bool outside = convex::outside_vertex_edge_voronoi_plane(p, b, a);
    BOOST_CHECK( outside );
  }
  // Back side of B voronoi plane
  {
    V p = V::make( 1.0, 1.0,  1.0);
    bool outside = convex::outside_vertex_edge_voronoi_plane(p, b, a);
    BOOST_CHECK( !outside );
  }
  // In B voronoi plane
  {
    V p = V::make( 0.0, 1.0,  1.0);
    bool outside = convex::outside_vertex_edge_voronoi_plane(p, b, a);
    BOOST_CHECK( outside );
  }

  // Second we use a test point that lies on the line

  // Front side of A voronoi plane
  {
    V p = V::make( 2.0, 0.0,  0.0);
    bool outside = convex::outside_vertex_edge_voronoi_plane(p, a, b);
    BOOST_CHECK( outside );
  }
  // Back side of A voronoi plane
  {
    V p = V::make( 0.0, 0.0,  0.0);
    bool outside = convex::outside_vertex_edge_voronoi_plane(p, a, b);
    BOOST_CHECK( !outside );
  }
  // In A voronoi plane
  {
    V p = V::make( 1.0, 0.0,  0.0);
    bool outside = convex::outside_vertex_edge_voronoi_plane(p, a, b);
    BOOST_CHECK( outside );
  }

  // Front side of B voronoi plane
  {
    V p = V::make( -1.0, 0.0,  0.0);
    bool outside = convex::outside_vertex_edge_voronoi_plane(p, b, a);
    BOOST_CHECK( outside );
  }
  // Back side of B voronoi plane
  {
    V p = V::make( 1.0, 0.0,  0.0);
    bool outside = convex::outside_vertex_edge_voronoi_plane(p, b, a);
    BOOST_CHECK( !outside );
  }
  // In B voronoi plane
  {
    V p = V::make( 0.0, 0.0,  0.0);
    bool outside = convex::outside_vertex_edge_voronoi_plane(p, b, a );
    BOOST_CHECK( outside );
  }
}

BOOST_AUTO_TEST_SUITE_END();
