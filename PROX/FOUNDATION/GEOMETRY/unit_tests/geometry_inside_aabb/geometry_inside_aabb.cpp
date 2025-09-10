#include <geometry.h>
#include <eigenhelperall.h>

#define BOOST_AUTO_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include <boost/test/unit_test_suite.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/test/test_tools.hpp>

BOOST_AUTO_TEST_SUITE(geometry);

BOOST_AUTO_TEST_CASE(inside_aabb_test)
{

    {
        using T = float;
        EigenVector3<T> const min_coord = EigenVector3<T>(-1.0, -1.0, -1.0);
        EigenVector3<T> const max_coord = EigenVector3<T>(1.0, 1.0, 1.0);

        geometry::AABBEigen<T> A = geometry::make_aabb(min_coord, max_coord);

        EigenVector3<T> const p0 = EigenVector3<T>(0.0, 0.0, 0.0);
        EigenVector3<T> const p1 = EigenVector3<T>(1.0, 1.0, 1.0);
        EigenVector3<T> const p2 = EigenVector3<T>(2.0, 2.0, 2.0);

        BOOST_CHECK(geometry::inside_aabb(p0, A) == true);
        BOOST_CHECK(geometry::inside_aabb(p1, A) == true);
        BOOST_CHECK(geometry::inside_aabb(p2, A) == false);
    }
}

BOOST_AUTO_TEST_SUITE_END();
