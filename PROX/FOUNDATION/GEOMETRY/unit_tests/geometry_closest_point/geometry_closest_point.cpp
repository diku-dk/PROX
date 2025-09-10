#include <geometry.h>
#include <eigenhelperall.h>

#define BOOST_AUTO_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include <boost/test/unit_test_suite.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/test/test_tools.hpp>

#include <vector>

BOOST_AUTO_TEST_SUITE(geometry);

BOOST_AUTO_TEST_CASE(closest_point_test)
{
    using T = float;

    {
        geometry::Line<T> const L = geometry::Line<T>(
            EigenVector3<T>(0, 0, 0), EigenVector3<T>(1, 0, 0)); // x-axis line

        EigenVector3<T> const q
            = geometry::closest_point_on_line(EigenVector3<T>(10, 4, 5), L);

        BOOST_CHECK_CLOSE(q(0), 10.0f, 0.01f);
        BOOST_CHECK_CLOSE(q(1), 0.0f, 0.01f);
        BOOST_CHECK_CLOSE(q(2), 0.0f, 0.01f);
    }
    {
        geometry::Plane<T> const P
            = geometry::make_plane<T>(EigenVector3<T>(1, 0, 0), 0); // y-z plane

        EigenVector3<T> const q
            = geometry::closest_point_on_plane(EigenVector3<T>(10, 4, 5), P);

        BOOST_CHECK_CLOSE(q(0), 0.0f, 0.01f);
        BOOST_CHECK_CLOSE(q(1), 4.0f, 0.01f);
        BOOST_CHECK_CLOSE(q(2), 5.0f, 0.01f);
    }
    {
        geometry::Line<T> const L1 = geometry::Line<T>(
            EigenVector3<T>(0, 0, 1), EigenVector3<T>(0, 1, 0)); // y-axis line
        geometry::Line<T> const L2 = geometry::Line<T>(
            EigenVector3<T>(0, 0, 0), EigenVector3<T>(1, 0, 0)); // x-axis line

        EigenVector3<T> p1, p2;

        geometry::closest_points_line_line_eigen(L1, L2, p1, p2);

        BOOST_CHECK_CLOSE(p1(0), 0.0f, 0.01f);
        BOOST_CHECK_CLOSE(p1(1), 0.0f, 0.01f);
        BOOST_CHECK_CLOSE(p1(2), 1.0f, 0.01f);

        BOOST_CHECK_CLOSE(p2(0), 0.0f, 0.01f);
        BOOST_CHECK_CLOSE(p2(1), 0.0f, 0.01f);
        BOOST_CHECK_CLOSE(p2(2), 0.0f, 0.01f);
    }
}

BOOST_AUTO_TEST_SUITE_END();
