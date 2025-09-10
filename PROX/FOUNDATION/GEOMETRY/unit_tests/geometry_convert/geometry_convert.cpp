#include <geometry.h>
#include <eigenhelperall.h>

#define BOOST_AUTO_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include <boost/test/unit_test_suite.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/test/test_tools.hpp>

#include <vector>

BOOST_AUTO_TEST_SUITE(geometry);

BOOST_AUTO_TEST_CASE(convert_test)
{
    using T = float;
    {
        EigenVector3<T> const center = EigenVector3<T>(0, 0, 0);
        T const radius = 1;
        geometry::Sphere<T> const sphere
            = geometry::make_sphere(center, radius);
        geometry::DOP<T, 6> const dop = geometry::convert<6, T>(sphere);

        BOOST_CHECK_EQUAL(dop.size(), 6u);
        BOOST_CHECK_CLOSE(dop(0).lower(), -1, 0.01);
        BOOST_CHECK_CLOSE(dop(1).lower(), -1, 0.01);
        BOOST_CHECK_CLOSE(dop(2).lower(), -1, 0.01);
        BOOST_CHECK_CLOSE(dop(0).upper(), 1, 0.01);
        BOOST_CHECK_CLOSE(dop(1).upper(), 1, 0.01);
        BOOST_CHECK_CLOSE(dop(2).upper(), 1, 0.01);
    }

    {
        EigenVector3<T> const min_coord = EigenVector3<T>(-1, -1, -1);
        EigenVector3<T> const max_coord = EigenVector3<T>(1, 1, 1);
        geometry::AABBEigen<T> const aabb
            = geometry::make_aabb(min_coord, max_coord);
        geometry::DOP<T, 6> const dop = geometry::convert<6, T>(aabb);

        BOOST_CHECK_EQUAL(dop.size(), 6u);
        BOOST_CHECK_CLOSE(dop(0).lower(), -1, 0.01);
        BOOST_CHECK_CLOSE(dop(1).lower(), -1, 0.01);
        BOOST_CHECK_CLOSE(dop(2).lower(), -1, 0.01);
        BOOST_CHECK_CLOSE(dop(0).upper(), 1, 0.01);
        BOOST_CHECK_CLOSE(dop(1).upper(), 1, 0.01);
        BOOST_CHECK_CLOSE(dop(2).upper(), 1, 0.01);
    }
}

BOOST_AUTO_TEST_SUITE_END();
