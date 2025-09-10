#include <geometry.h>
#include <tiny.h>

#define BOOST_AUTO_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include <boost/test/unit_test_suite.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/test/test_tools.hpp>

#include <vector>

BOOST_AUTO_TEST_SUITE(geometry);

BOOST_AUTO_TEST_CASE(line)
{
    using T = float;

    {
        EigenVector3<T> const p0 = EigenVector3<T>(0.0f, 0.0f, 0.0f);
        EigenVector3<T> const p1 = EigenVector3<T>(2.0f, 0.0f, 0.0f);

        geometry::Line<T> L = geometry::make_line(p0, p1);

        BOOST_CHECK_EQUAL(L.point()(0), 0.0f);
        BOOST_CHECK_EQUAL(L.point()(1), 0.0f);
        BOOST_CHECK_EQUAL(L.point()(2), 0.0f);

        BOOST_CHECK_EQUAL(L.direction()(0), 1.0f);
        BOOST_CHECK_EQUAL(L.direction()(1), 0.0f);
        BOOST_CHECK_EQUAL(L.direction()(2), 0.0f);
    }
    {
        EigenVector3<T> const p0 = EigenVector3<T>(0.0f, 0.0f, 0.0f);
        EigenVector3<T> const p1 = EigenVector3<T>(2.0f, 0.0f, 0.0f);

        geometry::Line<T> L
            = geometry::make_line(p0, p1, geometry::FROM_POINTS());

        BOOST_CHECK_EQUAL(L.point()(0), 0.0f);
        BOOST_CHECK_EQUAL(L.point()(1), 0.0f);
        BOOST_CHECK_EQUAL(L.point()(2), 0.0f);

        BOOST_CHECK_EQUAL(L.direction()(0), 1.0f);
        BOOST_CHECK_EQUAL(L.direction()(1), 0.0f);
        BOOST_CHECK_EQUAL(L.direction()(2), 0.0f);
    }
    {
        EigenVector3<T> const p = EigenVector3<T>(0.0f, 0.0f, 0.0f);
        EigenVector3<T> const d = EigenVector3<T>(0.5f, 0.0f, 0.0f);

        geometry::Line<T> L
            = geometry::make_line(p, d, geometry::FROM_DIRECTION());

        BOOST_CHECK_EQUAL(L.point()(0), 0.0f);
        BOOST_CHECK_EQUAL(L.point()(1), 0.0f);
        BOOST_CHECK_EQUAL(L.point()(2), 0.0f);

        BOOST_CHECK_EQUAL(L.direction()(0), 1.0f);
        BOOST_CHECK_EQUAL(L.direction()(1), 0.0f);
        BOOST_CHECK_EQUAL(L.direction()(2), 0.0f);
    }
}

BOOST_AUTO_TEST_SUITE_END();
