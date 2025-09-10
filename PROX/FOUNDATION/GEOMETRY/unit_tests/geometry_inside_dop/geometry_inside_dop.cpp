#include <geometry.h>
#include <tiny_math_types.h>

#define BOOST_AUTO_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include <boost/test/unit_test_suite.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/test/test_tools.hpp>

using MT = tiny::MathTypes<double>;
using V = MT::vector3_type;
using T = MT::real_type;
using VT = MT::value_traits;

BOOST_AUTO_TEST_SUITE(geometry);

BOOST_AUTO_TEST_CASE(inside_dop_test)
{
    std::vector<EigenVector3<T>> corners(2u);

    corners[0] = EigenVector3<T>(-1.0, -1.0, -1.0);
    corners[1] = EigenVector3<T>(1.0, 1.0, 1.0);

    geometry::DOP<T, 6u> const dop = geometry::make_dop(
        corners.begin(), corners.end(), geometry::make3<T>());

    BOOST_CHECK(geometry::is_valid(dop));

    {
        EigenVector3<T> const p = EigenVector3<T>(0.0, 0.0, 0.0);
        bool const test = geometry::inside_dop(p, dop);
        BOOST_CHECK(test);
    }
    {
        EigenVector3<T> const p = EigenVector3<T>(1.0, 1.0, 1.0);
        bool const test = geometry::inside_dop(p, dop);
        BOOST_CHECK(test);
    }
    {
        EigenVector3<T> const p = EigenVector3<T>(2.0, 0.0, 0.0);
        bool const test = geometry::inside_dop(p, dop);
        BOOST_CHECK(!test);
    }
}

BOOST_AUTO_TEST_CASE(outside_dop_test)
{
    std::vector<EigenVector3<T>> corners(2u);

    corners[0] = EigenVector3<T>(-1.0, -1.0, -1.0);
    corners[1] = EigenVector3<T>(1.0, 1.0, 1.0);

    geometry::DOP<T, 6u> const dop = geometry::make_dop(
        corners.begin(), corners.end(), geometry::make3<T>());

    BOOST_CHECK(geometry::is_valid(dop));

    {
        EigenVector3<T> const p = EigenVector3<T>(0.0, 0.0, 0.0);
        bool const test = geometry::outside_dop(p, dop);
        BOOST_CHECK(!test);
    }
    {
        EigenVector3<T> const p = EigenVector3<T>(1.0, 1.0, 1.0);
        bool const test = geometry::outside_dop(p, dop);
        BOOST_CHECK(!test);
    }
    {
        EigenVector3<T> const p = EigenVector3<T>(2.0, 0.0, 0.0);
        bool const test = geometry::outside_dop(p, dop);
        BOOST_CHECK(test);
    }
}

BOOST_AUTO_TEST_SUITE_END();
