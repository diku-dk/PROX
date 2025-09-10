#include <eigenhelperall.h>
#include <convex.h>

#define BOOST_AUTO_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include <boost/test/unit_test_suite.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/test/test_tools.hpp>

BOOST_AUTO_TEST_SUITE(convex_signed_distance_to_edge_face_vp);

BOOST_AUTO_TEST_CASE(case_by_case_test)
{

    using T = double;
    using real_type = T;

    EigenVector3<T> a = EigenVector3<T>(0.0, 0.0, 0.0);
    EigenVector3<T> b = EigenVector3<T>(1.0, 0.0, 0.0);
    ;
    EigenVector3<T> c = EigenVector3<T>(0.0, 1.0, 0.0);
    ;

  // Front side of AB voronoi plane
    {
        EigenVector3<T> p = EigenVector3<T>(-0.5, -1.0, 1.0);
        real_type sign_p = 0.0;
        sign_p = convex::signed_distance_to_edge_face_voronoi_plane(p, a, b, c);
        BOOST_CHECK_CLOSE(sign_p, 1.0, 0.01);

        sign_p = convex::signed_distance_to_edge_face_voronoi_plane(p, b, a, c);
        BOOST_CHECK_CLOSE(sign_p, 1.0, 0.01);
    }
  // Back side of AB voronoi plane
    {
        EigenVector3<T> p = EigenVector3<T>(-0.5, 1.0, 1.0);
        real_type sign_p = 0.0;
        sign_p = convex::signed_distance_to_edge_face_voronoi_plane(p, a, b, c);
        BOOST_CHECK_CLOSE(sign_p, -1.0, 0.01);

        sign_p = convex::signed_distance_to_edge_face_voronoi_plane(p, b, a, c);
        BOOST_CHECK_CLOSE(sign_p, -1.0, 0.01);
    }
  // In AB voronoi plane
    {
        EigenVector3<T> p = EigenVector3<T>(-0.5, 0.0, 1.0);
        real_type sign_p = 0.0;
        sign_p = convex::signed_distance_to_edge_face_voronoi_plane(p, a, b, c);
        BOOST_CHECK_CLOSE(sign_p, 0.0, 0.01);

        sign_p = convex::signed_distance_to_edge_face_voronoi_plane(p, b, a, c);
        BOOST_CHECK_CLOSE(sign_p, 0.0, 0.01);
    }
  // Front side of AC voronoi plane
    {
        EigenVector3<T> p = EigenVector3<T>(-1.0, 0.5, 1.0);
        real_type sign_p = 0.0;
        sign_p = convex::signed_distance_to_edge_face_voronoi_plane(p, a, c, b);
        BOOST_CHECK_CLOSE(sign_p, 1.0, 0.01);

        sign_p = convex::signed_distance_to_edge_face_voronoi_plane(p, c, a, b);
        BOOST_CHECK_CLOSE(sign_p, 1.0, 0.01);
    }
  // Back side of AC voronoi plane
    {
        EigenVector3<T> p = EigenVector3<T>(1.0, 0.5, 1.0);
        real_type sign_p = 0.0;
        sign_p = convex::signed_distance_to_edge_face_voronoi_plane(p, a, c, b);
        BOOST_CHECK_CLOSE(sign_p, -1.0, 0.01);

        sign_p = convex::signed_distance_to_edge_face_voronoi_plane(p, c, a, b);
        BOOST_CHECK_CLOSE(sign_p, -1.0, 0.01);
    }
  // In AC voronoi plane
    {
        EigenVector3<T> p = EigenVector3<T>(0.0, 0.5, 1.0);
        real_type sign_p = 0.0;
        sign_p = convex::signed_distance_to_edge_face_voronoi_plane(p, a, c, b);
        BOOST_CHECK_CLOSE(sign_p, 0.0, 0.01);

        sign_p = convex::signed_distance_to_edge_face_voronoi_plane(p, c, a, b);
        BOOST_CHECK_CLOSE(sign_p, 0.0, 0.01);
    }
  // Front side of BC voronoi plane
    {
        EigenVector3<T> p = EigenVector3<T>(1.0, 1.0, 1.0);
        real_type sign_p = 0.0;
        sign_p = convex::signed_distance_to_edge_face_voronoi_plane(p, b, c, a);
        BOOST_CHECK_CLOSE(sign_p, 0.70710678118654752440084436210485, 0.01);

        sign_p = convex::signed_distance_to_edge_face_voronoi_plane(p, c, b, a);
        BOOST_CHECK_CLOSE(sign_p, 0.70710678118654752440084436210485, 0.01);
    }
  // Back side of BC voronoi plane
    {
        EigenVector3<T> p = EigenVector3<T>(0.0, 0.0, 1.0);
        real_type sign_p = 0.0;
        sign_p = convex::signed_distance_to_edge_face_voronoi_plane(p, b, c, a);
        BOOST_CHECK_CLOSE(sign_p, -0.70710678118654752440084436210485, 0.01);

        sign_p = convex::signed_distance_to_edge_face_voronoi_plane(p, c, b, a);
        BOOST_CHECK_CLOSE(sign_p, -0.70710678118654752440084436210485, 0.01);
    }
  // In BC voronoi plane
    {
        EigenVector3<T> p = EigenVector3<T>(0.5, 0.5, 1.0);
        real_type sign_p = 0.0;
        sign_p = convex::signed_distance_to_edge_face_voronoi_plane(p, b, c, a);
        BOOST_CHECK_CLOSE(sign_p, 0.0, 0.01);

        sign_p = convex::signed_distance_to_edge_face_voronoi_plane(p, c, b, a);
        BOOST_CHECK_CLOSE(sign_p, 0.0, 0.01);
    }

  // We just used a point above the face plane, next we will try using a test point lying in the face plane

  // Front side of AB voronoi plane
    {
        EigenVector3<T> p = EigenVector3<T>(-0.5, -1.0, 0.0);
        real_type sign_p = 0.0;
        sign_p = convex::signed_distance_to_edge_face_voronoi_plane(p, a, b, c);
        BOOST_CHECK_CLOSE(sign_p, 1.0, 0.01);

        sign_p = convex::signed_distance_to_edge_face_voronoi_plane(p, b, a, c);
        BOOST_CHECK_CLOSE(sign_p, 1.0, 0.01);
    }
  // Back side of AB voronoi plane
    {
        EigenVector3<T> p = EigenVector3<T>(-0.5, 1.0, 0.0);
        real_type sign_p = 0.0;
        sign_p = convex::signed_distance_to_edge_face_voronoi_plane(p, a, b, c);
        BOOST_CHECK_CLOSE(sign_p, -1.0, 0.01);

        sign_p = convex::signed_distance_to_edge_face_voronoi_plane(p, b, a, c);
        BOOST_CHECK_CLOSE(sign_p, -1.0, 0.01);
    }
  // In AB voronoi plane
    {
        EigenVector3<T> p = EigenVector3<T>(-0.5, 0.0, 0.0);
        real_type sign_p = 0.0;
        sign_p = convex::signed_distance_to_edge_face_voronoi_plane(p, a, b, c);
        BOOST_CHECK_CLOSE(sign_p, 0.0, 0.01);

        sign_p = convex::signed_distance_to_edge_face_voronoi_plane(p, b, a, c);
        BOOST_CHECK_CLOSE(sign_p, 0.0, 0.01);
    }
  // Front side of AC voronoi plane
    {
        EigenVector3<T> p = EigenVector3<T>(-1.0, 0.5, 0.0);
        real_type sign_p = 0.0;
        sign_p = convex::signed_distance_to_edge_face_voronoi_plane(p, a, c, b);
        BOOST_CHECK_CLOSE(sign_p, 1.0, 0.01);

        sign_p = convex::signed_distance_to_edge_face_voronoi_plane(p, c, a, b);
        BOOST_CHECK_CLOSE(sign_p, 1.0, 0.01);
    }
  // Back side of AC voronoi plane
    {
        EigenVector3<T> p = EigenVector3<T>(1.0, 0.5, 0.0);
        real_type sign_p = 0.0;
        sign_p = convex::signed_distance_to_edge_face_voronoi_plane(p, a, c, b);
        BOOST_CHECK_CLOSE(sign_p, -1.0, 0.01);

        sign_p = convex::signed_distance_to_edge_face_voronoi_plane(p, c, a, b);
        BOOST_CHECK_CLOSE(sign_p, -1.0, 0.01);
    }
  // In AC voronoi plane
    {
        EigenVector3<T> p = EigenVector3<T>(0.0, 0.5, 0.0);
        real_type sign_p = 0.0;
        sign_p = convex::signed_distance_to_edge_face_voronoi_plane(p, a, c, b);
        BOOST_CHECK_CLOSE(sign_p, 0.0, 0.01);

        sign_p = convex::signed_distance_to_edge_face_voronoi_plane(p, c, a, b);
        BOOST_CHECK_CLOSE(sign_p, 0.0, 0.01);
    }
  // Front side of BC voronoi plane
    {
        EigenVector3<T> p = EigenVector3<T>(1.0, 1.0, 0.0);
        real_type sign_p = 0.0;
        sign_p = convex::signed_distance_to_edge_face_voronoi_plane(p, b, c, a);
        BOOST_CHECK_CLOSE(sign_p, 0.70710678118654752440084436210485, 0.01);

        sign_p = convex::signed_distance_to_edge_face_voronoi_plane(p, c, b, a);
        BOOST_CHECK_CLOSE(sign_p, 0.70710678118654752440084436210485, 0.01);
    }
  // Back side of BC voronoi plane
    {
        EigenVector3<T> p = EigenVector3<T>(0.0, 0.0, 0.0);
        real_type sign_p = 0.0;
        sign_p = convex::signed_distance_to_edge_face_voronoi_plane(p, b, c, a);
        BOOST_CHECK_CLOSE(sign_p, -0.70710678118654752440084436210485, 0.01);

        sign_p = convex::signed_distance_to_edge_face_voronoi_plane(p, c, b, a);
        BOOST_CHECK_CLOSE(sign_p, -0.70710678118654752440084436210485, 0.01);
    }
  // In BC voronoi plane
    {
        EigenVector3<T> p = EigenVector3<T>(0.5, 0.5, 0.0);
        real_type sign_p = 0.0;
        sign_p = convex::signed_distance_to_edge_face_voronoi_plane(p, b, c, a);
        BOOST_CHECK_CLOSE(sign_p, 0.0, 0.01);

        sign_p = convex::signed_distance_to_edge_face_voronoi_plane(p, c, b, a);
        BOOST_CHECK_CLOSE(sign_p, 0.0, 0.01);
    }

  // Finally we will use a test point lying below the face-plane

  // Front side of AB voronoi plane
    {
        EigenVector3<T> p = EigenVector3<T>(-0.5, -1.0, -1.0);
        real_type sign_p = 0.0;
        sign_p = convex::signed_distance_to_edge_face_voronoi_plane(p, a, b, c);
        BOOST_CHECK_CLOSE(sign_p, 1.0, 0.01);

        sign_p = convex::signed_distance_to_edge_face_voronoi_plane(p, b, a, c);
        BOOST_CHECK_CLOSE(sign_p, 1.0, 0.01);
    }
  // Back side of AB voronoi plane
    {
        EigenVector3<T> p = EigenVector3<T>(-0.5, 1.0, -1.0);
        real_type sign_p = 0.0;
        sign_p = convex::signed_distance_to_edge_face_voronoi_plane(p, a, b, c);
        BOOST_CHECK_CLOSE(sign_p, -1.0, 0.01);

        sign_p = convex::signed_distance_to_edge_face_voronoi_plane(p, b, a, c);
        BOOST_CHECK_CLOSE(sign_p, -1.0, 0.01);
    }
  // In AB voronoi plane
    {
        EigenVector3<T> p = EigenVector3<T>(-0.5, 0.0, -1.0);
        real_type sign_p = 0.0;
        sign_p = convex::signed_distance_to_edge_face_voronoi_plane(p, a, b, c);
        BOOST_CHECK_CLOSE(sign_p, 0.0, 0.01);

        sign_p = convex::signed_distance_to_edge_face_voronoi_plane(p, b, a, c);
        BOOST_CHECK_CLOSE(sign_p, 0.0, 0.01);
    }
  // Front side of AC voronoi plane
    {
        EigenVector3<T> p = EigenVector3<T>(-1.0, 0.5, -1.0);
        real_type sign_p = 0.0;
        sign_p = convex::signed_distance_to_edge_face_voronoi_plane(p, a, c, b);
        BOOST_CHECK_CLOSE(sign_p, 1.0, 0.01);

        sign_p = convex::signed_distance_to_edge_face_voronoi_plane(p, c, a, b);
        BOOST_CHECK_CLOSE(sign_p, 1.0, 0.01);
    }
  // Back side of AC voronoi plane
    {
        EigenVector3<T> p = EigenVector3<T>(1.0, 0.5, -1.0);
        real_type sign_p = 0.0;
        sign_p = convex::signed_distance_to_edge_face_voronoi_plane(p, a, c, b);
        BOOST_CHECK_CLOSE(sign_p, -1.0, 0.01);

        sign_p = convex::signed_distance_to_edge_face_voronoi_plane(p, c, a, b);
        BOOST_CHECK_CLOSE(sign_p, -1.0, 0.01);
    }
  // In AC voronoi plane
    {
        EigenVector3<T> p = EigenVector3<T>(0.0, 0.5, -1.0);
        real_type sign_p = 0.0;
        sign_p = convex::signed_distance_to_edge_face_voronoi_plane(p, a, c, b);
        BOOST_CHECK_CLOSE(sign_p, 0.0, 0.01);

        sign_p = convex::signed_distance_to_edge_face_voronoi_plane(p, c, a, b);
        BOOST_CHECK_CLOSE(sign_p, 0.0, 0.01);
    }
  // Front side of BC voronoi plane
    {
        EigenVector3<T> p = EigenVector3<T>(1.0, 1.0, -1.0);
        real_type sign_p = 0.0;
        sign_p = convex::signed_distance_to_edge_face_voronoi_plane(p, b, c, a);
        BOOST_CHECK_CLOSE(sign_p, 0.70710678118654752440084436210485, 0.01);

        sign_p = convex::signed_distance_to_edge_face_voronoi_plane(p, c, b, a);
        BOOST_CHECK_CLOSE(sign_p, 0.70710678118654752440084436210485, 0.01);
    }
  // Back side of BC voronoi plane
    {
        EigenVector3<T> p = EigenVector3<T>(0.0, 0.0, -1.0);
        real_type sign_p = 0.0;
        sign_p = convex::signed_distance_to_edge_face_voronoi_plane(p, b, c, a);
        BOOST_CHECK_CLOSE(sign_p, -0.70710678118654752440084436210485, 0.01);

        sign_p = convex::signed_distance_to_edge_face_voronoi_plane(p, c, b, a);
        BOOST_CHECK_CLOSE(sign_p, -0.70710678118654752440084436210485, 0.01);
    }
  // In BC voronoi plane
    {
        EigenVector3<T> p = EigenVector3<T>(0.5, 0.5, -1.0);
        real_type sign_p = 0.0;
        sign_p = convex::signed_distance_to_edge_face_voronoi_plane(p, b, c, a);
        BOOST_CHECK_CLOSE(sign_p, 0.0, 0.01);

        sign_p = convex::signed_distance_to_edge_face_voronoi_plane(p, c, b, a);
        BOOST_CHECK_CLOSE(sign_p, 0.0, 0.01);
    }
}

BOOST_AUTO_TEST_SUITE_END();
