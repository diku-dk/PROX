#include <geometry.h>
#include <tiny.h>

#define BOOST_AUTO_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include <boost/test/unit_test_suite.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/test/test_tools.hpp>

#include <vector>

BOOST_AUTO_TEST_SUITE(geometry);

BOOST_AUTO_TEST_CASE(intersection)
{
    using T = float;

    {
        geometry::Plane<T> const A
            = geometry::make_plane<T>(EigenVector3<T>(0, 0, 1), -1.0);
        geometry::Plane<T> const B
            = geometry::make_plane<T>(EigenVector3<T>(0, 1, 0), -1.0);
        geometry::Line<T> const L = geometry::make_intersection(A, B);

        BOOST_CHECK_CLOSE(geometry::get_distance(L.point(), A), 0.0, 0.01f);
        BOOST_CHECK_CLOSE(geometry::get_distance(L.point(), B), 0.0, 0.01f);
        BOOST_CHECK_CLOSE(dot(L.direction(), A.normal()), 0.0, 0.01f);
        BOOST_CHECK_CLOSE(dot(L.direction(), B.normal()), 0.0, 0.01f);
        BOOST_CHECK_CLOSE(dot(L.direction(), L.direction()), 1.0, 0.01f);
    }

    {
        geometry::Plane<T> const A
            = geometry::make_plane<T>(EigenVector3<T>(0, 0, 1), -1.0);
        geometry::Plane<T> const B
            = geometry::make_plane<T>(EigenVector3<T>(0, 1, 0), -1.0);
        geometry::Plane<T> const C
            = geometry::make_plane<T>(EigenVector3<T>(1, 0, 0), -1.0);

        const EigenVector3<T> p = geometry::make_intersection(A, B, C);

        BOOST_CHECK_CLOSE(geometry::get_distance(p, A), 0.0, 0.01f);
        BOOST_CHECK_CLOSE(geometry::get_distance(p, B), 0.0, 0.01f);
        BOOST_CHECK_CLOSE(geometry::get_distance(p, C), 0.0, 0.01f);
    }
}

BOOST_AUTO_TEST_SUITE_END();
