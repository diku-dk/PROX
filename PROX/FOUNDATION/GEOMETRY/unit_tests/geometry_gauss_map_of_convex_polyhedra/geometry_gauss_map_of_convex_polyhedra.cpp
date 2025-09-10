#include <geometry.h>
#include <tiny.h>

#define BOOST_AUTO_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include <boost/test/unit_test_suite.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/test/test_tools.hpp>

using T = float;

BOOST_AUTO_TEST_SUITE(geometry);

BOOST_AUTO_TEST_CASE(gauss_map_tetrahedron_test)
{

    {
        geometry::TetrahedronEigen<T> A;

        geometry::GaussMapOfConvexPolyhedra<T> G = geometry::make_gauss_map(A);

    // Search for face features
        {
            EigenVector3<T> const s = unit(EigenVector3<T>(0.0, 0.0, -1.0));

            EigenVector3<T> n;
            std::vector<EigenVector3<T>> points;

            bool found = G.search_for_feature(s, points, n);

            BOOST_CHECK(found);
            BOOST_CHECK_EQUAL(points.size(), 3);

            BOOST_CHECK_CLOSE(points[0][0], 0.0, 0.01);
            BOOST_CHECK_CLOSE(points[0][1], 0.0, 0.01);
            BOOST_CHECK_CLOSE(points[0][2], 0.0, 0.01);

            BOOST_CHECK_CLOSE(points[1][0], 0.0, 0.01);
            BOOST_CHECK_CLOSE(points[1][1], 1.0, 0.01);
            BOOST_CHECK_CLOSE(points[1][2], 0.0, 0.01);

            BOOST_CHECK_CLOSE(points[2][0], 1.0, 0.01);
            BOOST_CHECK_CLOSE(points[2][1], 0.0, 0.01);
            BOOST_CHECK_CLOSE(points[2][2], 0.0, 0.01);

            BOOST_CHECK_CLOSE(n[0], 0.0, 0.01);
            BOOST_CHECK_CLOSE(n[1], 0.0, 0.01);
            BOOST_CHECK_CLOSE(n[2], -1.0, 0.01);
        }
        {
            EigenVector3<T> const s = unit(EigenVector3<T>(-1.0, 0.0, 0.0));

            EigenVector3<T> n;
            std::vector<EigenVector3<T>> points;

            bool found = G.search_for_feature(s, points, n);

            BOOST_CHECK(found);
            BOOST_CHECK_EQUAL(points.size(), 3);

            BOOST_CHECK_CLOSE(points[0][0], 0.0, 0.01);
            BOOST_CHECK_CLOSE(points[0][1], 1.0, 0.01);
            BOOST_CHECK_CLOSE(points[0][2], 0.0, 0.01);

            BOOST_CHECK_CLOSE(points[1][0], 0.0, 0.01);
            BOOST_CHECK_CLOSE(points[1][1], 0.0, 0.01);
            BOOST_CHECK_CLOSE(points[1][2], 0.0, 0.01);

            BOOST_CHECK_CLOSE(points[2][0], 0.0, 0.01);
            BOOST_CHECK_CLOSE(points[2][1], 0.0, 0.01);
            BOOST_CHECK_CLOSE(points[2][2], 1.0, 0.01);

            BOOST_CHECK_CLOSE(n[0], -1.0, 0.01);
            BOOST_CHECK_CLOSE(n[1], 0.0, 0.01);
            BOOST_CHECK_CLOSE(n[2], 0.0, 0.01);
        }
        {
            EigenVector3<T> const s = unit(EigenVector3<T>(0.0, -1.0, 0.0));

            EigenVector3<T> n;
            std::vector<EigenVector3<T>> points;

            bool found = G.search_for_feature(s, points, n);

            BOOST_CHECK(found);
            BOOST_CHECK_EQUAL(points.size(), 3);

            BOOST_CHECK_CLOSE(points[0][0], 0.0, 0.01);
            BOOST_CHECK_CLOSE(points[0][1], 0.0, 0.01);
            BOOST_CHECK_CLOSE(points[0][2], 0.0, 0.01);

            BOOST_CHECK_CLOSE(points[1][0], 1.0, 0.01);
            BOOST_CHECK_CLOSE(points[1][1], 0.0, 0.01);
            BOOST_CHECK_CLOSE(points[1][2], 0.0, 0.01);

            BOOST_CHECK_CLOSE(points[2][0], 0.0, 0.01);
            BOOST_CHECK_CLOSE(points[2][1], 0.0, 0.01);
            BOOST_CHECK_CLOSE(points[2][2], 1.0, 0.01);

            BOOST_CHECK_CLOSE(n[0], 0.0, 0.01);
            BOOST_CHECK_CLOSE(n[1], -1.0, 0.01);
            BOOST_CHECK_CLOSE(n[2], 0.0, 0.01);
        }
        {
            EigenVector3<T> const s = unit(EigenVector3<T>(1.0, 1.0, 1.0));

            EigenVector3<T> n;
            std::vector<EigenVector3<T>> points;

            bool found = G.search_for_feature(s, points, n);

            BOOST_CHECK(found);
            BOOST_CHECK_EQUAL(points.size(), 3);

            BOOST_CHECK_CLOSE(points[0][0], 1.0, 0.01);
            BOOST_CHECK_CLOSE(points[0][1], 0.0, 0.01);
            BOOST_CHECK_CLOSE(points[0][2], 0.0, 0.01);

            BOOST_CHECK_CLOSE(points[1][0], 0.0, 0.01);
            BOOST_CHECK_CLOSE(points[1][1], 1.0, 0.01);
            BOOST_CHECK_CLOSE(points[1][2], 0.0, 0.01);

            BOOST_CHECK_CLOSE(points[2][0], 0.0, 0.01);
            BOOST_CHECK_CLOSE(points[2][1], 0.0, 0.01);
            BOOST_CHECK_CLOSE(points[2][2], 1.0, 0.01);

            BOOST_CHECK_CLOSE(n[0], s[0], 0.01);
            BOOST_CHECK_CLOSE(n[1], s[1], 0.01);
            BOOST_CHECK_CLOSE(n[2], s[2], 0.01);
        }

    // Search for vertex features
        {
            EigenVector3<T> const s = unit(EigenVector3<T>(0.0, 0.0, 1.0));

            EigenVector3<T> n;
            std::vector<EigenVector3<T>> points;

            bool found = G.search_for_feature(s, points, n);

            BOOST_CHECK(found);
            BOOST_CHECK_EQUAL(points.size(), 1);

            BOOST_CHECK_CLOSE(points[0][0], 0.0, 0.01);
            BOOST_CHECK_CLOSE(points[0][1], 0.0, 0.01);
            BOOST_CHECK_CLOSE(points[0][2], 1.0, 0.01);

            BOOST_CHECK_CLOSE(n[0], s[0], 0.01);
            BOOST_CHECK_CLOSE(n[1], s[1], 0.01);
            BOOST_CHECK_CLOSE(n[2], s[2], 0.01);
        }
        {
            EigenVector3<T> const s = unit(EigenVector3<T>(0.0, 1.0, 0.0));

            EigenVector3<T> n;
            std::vector<EigenVector3<T>> points;

            bool found = G.search_for_feature(s, points, n);

            BOOST_CHECK(found);
            BOOST_CHECK_EQUAL(points.size(), 1);

            BOOST_CHECK_CLOSE(points[0][0], 0.0, 0.01);
            BOOST_CHECK_CLOSE(points[0][1], 1.0, 0.01);
            BOOST_CHECK_CLOSE(points[0][2], 0.0, 0.01);

            BOOST_CHECK_CLOSE(n[0], s[0], 0.01);
            BOOST_CHECK_CLOSE(n[1], s[1], 0.01);
            BOOST_CHECK_CLOSE(n[2], s[2], 0.01);
        }
        {
            EigenVector3<T> const s = unit(EigenVector3<T>(1.0, 0.0, 0.0));

            EigenVector3<T> n;
            std::vector<EigenVector3<T>> points;

            bool found = G.search_for_feature(s, points, n);

            BOOST_CHECK(found);
            BOOST_CHECK_EQUAL(points.size(), 1);

            BOOST_CHECK_CLOSE(points[0][0], 1.0, 0.01);
            BOOST_CHECK_CLOSE(points[0][1], 0.0, 0.01);
            BOOST_CHECK_CLOSE(points[0][2], 0.0, 0.01);

            BOOST_CHECK_CLOSE(n[0], s[0], 0.01);
            BOOST_CHECK_CLOSE(n[1], s[1], 0.01);
            BOOST_CHECK_CLOSE(n[2], s[2], 0.01);
        }
        {
            EigenVector3<T> const s = unit(EigenVector3<T>(-1.0, -1.0, -1.0));

            EigenVector3<T> n;
            std::vector<EigenVector3<T>> points;

            bool found = G.search_for_feature(s, points, n);

            BOOST_CHECK(found);
            BOOST_CHECK_EQUAL(points.size(), 1);

            BOOST_CHECK_CLOSE(points[0][0], 0.0, 0.01);
            BOOST_CHECK_CLOSE(points[0][1], 0.0, 0.01);
            BOOST_CHECK_CLOSE(points[0][2], 0.0, 0.01);

            BOOST_CHECK_CLOSE(n[0], s[0], 0.01);
            BOOST_CHECK_CLOSE(n[1], s[1], 0.01);
            BOOST_CHECK_CLOSE(n[2], s[2], 0.01);
        }

    // Search for edge features
        {
            EigenVector3<T> const s = unit(EigenVector3<T>(-1.0, -1.0, 0.0));

            EigenVector3<T> n;
            std::vector<EigenVector3<T>> points;

            bool found = G.search_for_feature(s, points, n);

            BOOST_CHECK(found);
            BOOST_CHECK_EQUAL(points.size(), 2);

            BOOST_CHECK_CLOSE(points[0][0], 0.0, 0.01);
            BOOST_CHECK_CLOSE(points[0][1], 0.0, 0.01);
            BOOST_CHECK_CLOSE(points[0][2], 0.0, 0.01);

            BOOST_CHECK_CLOSE(points[1][0], 0.0, 0.01);
            BOOST_CHECK_CLOSE(points[1][1], 0.0, 0.01);
            BOOST_CHECK_CLOSE(points[1][2], 1.0, 0.01);

            BOOST_CHECK_CLOSE(n[0], s[0], 0.01);
            BOOST_CHECK_CLOSE(n[1], s[1], 0.01);
            BOOST_CHECK_CLOSE(n[2], s[2], 0.01);
        }
        {
            EigenVector3<T> const s = unit(EigenVector3<T>(-1.0, 0.0, -1.0));

            EigenVector3<T> n;
            std::vector<EigenVector3<T>> points;

            bool found = G.search_for_feature(s, points, n);

            BOOST_CHECK(found);
            BOOST_CHECK_EQUAL(points.size(), 2);

            BOOST_CHECK_CLOSE(points[0][0], 0.0, 0.01);
            BOOST_CHECK_CLOSE(points[0][1], 0.0, 0.01);
            BOOST_CHECK_CLOSE(points[0][2], 0.0, 0.01);

            BOOST_CHECK_CLOSE(points[1][0], 0.0, 0.01);
            BOOST_CHECK_CLOSE(points[1][1], 1.0, 0.01);
            BOOST_CHECK_CLOSE(points[1][2], 0.0, 0.01);

            BOOST_CHECK_CLOSE(n[0], s[0], 0.01);
            BOOST_CHECK_CLOSE(n[1], s[1], 0.01);
            BOOST_CHECK_CLOSE(n[2], s[2], 0.01);
        }
        {
            EigenVector3<T> const s = unit(EigenVector3<T>(0.0, -1.0, -1.0));

            EigenVector3<T> n;
            std::vector<EigenVector3<T>> points;

            bool found = G.search_for_feature(s, points, n);

            BOOST_CHECK(found);
            BOOST_CHECK_EQUAL(points.size(), 2);

            BOOST_CHECK_CLOSE(points[0][0], 0.0, 0.01);
            BOOST_CHECK_CLOSE(points[0][1], 0.0, 0.01);
            BOOST_CHECK_CLOSE(points[0][2], 0.0, 0.01);

            BOOST_CHECK_CLOSE(points[1][0], 1.0, 0.01);
            BOOST_CHECK_CLOSE(points[1][1], 0.0, 0.01);
            BOOST_CHECK_CLOSE(points[1][2], 0.0, 0.01);

            BOOST_CHECK_CLOSE(n[0], s[0], 0.01);
            BOOST_CHECK_CLOSE(n[1], s[1], 0.01);
            BOOST_CHECK_CLOSE(n[2], s[2], 0.01);
        }
        {
            EigenVector3<T> const s = unit(EigenVector3<T>(1.0, 1.0, 0.5));

            EigenVector3<T> n;
            std::vector<EigenVector3<T>> points;

            bool found = G.search_for_feature(s, points, n);

            BOOST_CHECK(found);
            BOOST_CHECK_EQUAL(points.size(), 2);

            BOOST_CHECK_CLOSE(points[0][0], 1.0, 0.01);
            BOOST_CHECK_CLOSE(points[0][1], 0.0, 0.01);
            BOOST_CHECK_CLOSE(points[0][2], 0.0, 0.01);

            BOOST_CHECK_CLOSE(points[1][0], 0.0, 0.01);
            BOOST_CHECK_CLOSE(points[1][1], 1.0, 0.01);
            BOOST_CHECK_CLOSE(points[1][2], 0.0, 0.01);

            BOOST_CHECK_CLOSE(n[0], s[0], 0.01);
            BOOST_CHECK_CLOSE(n[1], s[1], 0.01);
            BOOST_CHECK_CLOSE(n[2], s[2], 0.01);
        }
        {
            EigenVector3<T> const s = unit(EigenVector3<T>(1.0, -0.5, 1.0));

            EigenVector3<T> n;
            std::vector<EigenVector3<T>> points;

            bool found = G.search_for_feature(s, points, n);

            BOOST_CHECK(found);
            BOOST_CHECK_EQUAL(points.size(), 2);

            BOOST_CHECK_CLOSE(points[0][0], 1.0, 0.01);
            BOOST_CHECK_CLOSE(points[0][1], 0.0, 0.01);
            BOOST_CHECK_CLOSE(points[0][2], 0.0, 0.01);

            BOOST_CHECK_CLOSE(points[1][0], 0.0, 0.01);
            BOOST_CHECK_CLOSE(points[1][1], 0.0, 0.01);
            BOOST_CHECK_CLOSE(points[1][2], 1.0, 0.01);

            BOOST_CHECK_CLOSE(n[0], s[0], 0.01);
            BOOST_CHECK_CLOSE(n[1], s[1], 0.01);
            BOOST_CHECK_CLOSE(n[2], s[2], 0.01);
        }
        {
            EigenVector3<T> const s = unit(EigenVector3<T>(-0.5, 1.0, 1.0));

            EigenVector3<T> n;
            std::vector<EigenVector3<T>> points;

            bool found = G.search_for_feature(s, points, n);

            BOOST_CHECK(found);
            BOOST_CHECK_EQUAL(points.size(), 2);

            BOOST_CHECK_CLOSE(points[0][0], 0.0, 0.01);
            BOOST_CHECK_CLOSE(points[0][1], 1.0, 0.01);
            BOOST_CHECK_CLOSE(points[0][2], 0.0, 0.01);

            BOOST_CHECK_CLOSE(points[1][0], 0.0, 0.01);
            BOOST_CHECK_CLOSE(points[1][1], 0.0, 0.01);
            BOOST_CHECK_CLOSE(points[1][2], 1.0, 0.01);

            BOOST_CHECK_CLOSE(n[0], s[0], 0.01);
            BOOST_CHECK_CLOSE(n[1], s[1], 0.01);
            BOOST_CHECK_CLOSE(n[2], s[2], 0.01);
        }
    }
}

BOOST_AUTO_TEST_SUITE_END();
