#include <tiny.h>
#include <convex.h>

#define BOOST_AUTO_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include <boost/test/unit_test_suite.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/test/test_tools.hpp>

BOOST_AUTO_TEST_SUITE(convex_signed_distance_to_vertex_edge_vp);

BOOST_AUTO_TEST_CASE(case_by_case_test)
{
    using T = double;
    using real_type = T;
    EigenVector3<T> a = EigenVector3<T>(1.0, 0.0, 0.0);
    EigenVector3<T> b = EigenVector3<T>(0.0, 0.0, 0.0);

  // First we use a test point that does not lie on the line

  // Front side of A voronoi plane
    {
        EigenVector3<T> p = EigenVector3<T>(2.0, 1.0, 1.0);
        real_type sign_p = 0.0;
        sign_p = convex::signed_distance_to_vertex_edge_voronoi_plane(p, a, b);
        BOOST_CHECK_CLOSE(sign_p, 1.0, 0.01);
    }
  // Back side of A voronoi plane
    {
        EigenVector3<T> p = EigenVector3<T>(0.0, 1.0, 1.0);
        real_type sign_p = 0.0;
        sign_p = convex::signed_distance_to_vertex_edge_voronoi_plane(p, a, b);
        BOOST_CHECK_CLOSE(sign_p, -1.0, 0.01);
    }
  // In A voronoi plane
    {
        EigenVector3<T> p = EigenVector3<T>(1.0, 1.0, 1.0);
        real_type sign_p = 0.0;
        sign_p = convex::signed_distance_to_vertex_edge_voronoi_plane(p, a, b);
        BOOST_CHECK_CLOSE(sign_p, 0.0, 0.01);
    }

  // Front side of B voronoi plane
    {
        EigenVector3<T> p = EigenVector3<T>(-1.0, 1.0, 1.0);
        real_type sign_p = 0.0;
        sign_p = convex::signed_distance_to_vertex_edge_voronoi_plane(p, b, a);
        BOOST_CHECK_CLOSE(sign_p, 1.0, 0.01);
    }
  // Back side of B voronoi plane
    {
        EigenVector3<T> p = EigenVector3<T>(1.0, 1.0, 1.0);
        real_type sign_p = 0.0;
        sign_p = convex::signed_distance_to_vertex_edge_voronoi_plane(p, b, a);
        BOOST_CHECK_CLOSE(sign_p, -1.0, 0.01);
    }
  // In B voronoi plane
    {
        EigenVector3<T> p = EigenVector3<T>(0.0, 1.0, 1.0);
        real_type sign_p = 0.0;
        sign_p = convex::signed_distance_to_vertex_edge_voronoi_plane(p, b, a);
        BOOST_CHECK_CLOSE(sign_p, 0.0, 0.01);
    }

  // Second we use a test point that lies on the line

  // Front side of A voronoi plane
    {
        EigenVector3<T> p = EigenVector3<T>(2.0, 0.0, 0.0);
        real_type sign_p = 0.0;
        sign_p = convex::signed_distance_to_vertex_edge_voronoi_plane(p, a, b);
        BOOST_CHECK_CLOSE(sign_p, 1.0, 0.01);
    }
  // Back side of A voronoi plane
    {
        EigenVector3<T> p = EigenVector3<T>(0.0, 0.0, 0.0);
        real_type sign_p = 0.0;
        sign_p = convex::signed_distance_to_vertex_edge_voronoi_plane(p, a, b);
        BOOST_CHECK_CLOSE(sign_p, -1.0, 0.01);
    }
  // In A voronoi plane
    {
        EigenVector3<T> p = EigenVector3<T>(1.0, 0.0, 0.0);
        real_type sign_p = 0.0;
        sign_p = convex::signed_distance_to_vertex_edge_voronoi_plane(p, a, b);
        BOOST_CHECK_CLOSE(sign_p, 0.0, 0.01);
    }

  // Front side of B voronoi plane
    {
        EigenVector3<T> p = EigenVector3<T>(-1.0, 0.0, 0.0);
        real_type sign_p = 0.0;
        sign_p = convex::signed_distance_to_vertex_edge_voronoi_plane(p, b, a);
        BOOST_CHECK_CLOSE(sign_p, 1.0, 0.01);
    }
  // Back side of B voronoi plane
    {
        EigenVector3<T> p = EigenVector3<T>(1.0, 0.0, 0.0);
        real_type sign_p = 0.0;
        sign_p = convex::signed_distance_to_vertex_edge_voronoi_plane(p, b, a);
        BOOST_CHECK_CLOSE(sign_p, -1.0, 0.01);
    }
  // In B voronoi plane
    {
        EigenVector3<T> p = EigenVector3<T>(0.0, 0.0, 0.0);
        real_type sign_p = 0.0;
        sign_p = convex::signed_distance_to_vertex_edge_voronoi_plane(p, b, a);
        BOOST_CHECK_CLOSE(sign_p, 0.0, 0.01);
    }
}

BOOST_AUTO_TEST_SUITE_END();
