#include <geometry.h>
#include <eigenhelperall.h>

#define BOOST_AUTO_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include <boost/test/unit_test_suite.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/test/test_tools.hpp>

using T = double;

BOOST_AUTO_TEST_SUITE(geometry);

BOOST_AUTO_TEST_CASE(aabb_test)
{

    {
        geometry::AABBEigen<T> A;
        BOOST_CHECK(geometry::is_valid(A) == true);

        BOOST_CHECK_CLOSE(A.min()(0), -0.5, 0.01);
        BOOST_CHECK_CLOSE(A.min()(1), -0.5, 0.01);
        BOOST_CHECK_CLOSE(A.min()(2), -0.5, 0.01);

        BOOST_CHECK_CLOSE(A.max()(0), 0.5, 0.01);
        BOOST_CHECK_CLOSE(A.max()(1), 0.5, 0.01);
        BOOST_CHECK_CLOSE(A.max()(2), 0.5, 0.01);
    }

    {
        EigenVector3<T> const min_coord = EigenVector3<T>(-1.0, -1.0, -1.0);
        EigenVector3<T> const max_coord = EigenVector3<T>(1.0, 1.0, 1.0);

        geometry::AABBEigen<T> A = geometry::make_aabb(min_coord, max_coord);

        BOOST_CHECK(geometry::is_valid(A) == true);

        BOOST_CHECK_CLOSE(A.min()(0), -1.0, 0.01);
        BOOST_CHECK_CLOSE(A.min()(1), -1.0, 0.01);
        BOOST_CHECK_CLOSE(A.min()(2), -1.0, 0.01);

        BOOST_CHECK_CLOSE(A.max()(0), 1.0, 0.01);
        BOOST_CHECK_CLOSE(A.max()(1), 1.0, 0.01);
        BOOST_CHECK_CLOSE(A.max()(2), 1.0, 0.01);
    }

    {
        EigenVector3<T> const min_coord = EigenVector3<T>(-1.0, -1.0, -1.0);
        EigenVector3<T> const max_coord = EigenVector3<T>(1.0, 1.0, 1.0);

        geometry::AABBEigen<T> A;

        A.min() = min_coord;
        A.max() = max_coord;

        BOOST_CHECK(geometry::is_valid(A) == true);

        BOOST_CHECK_CLOSE(A.min()(0), -1.0, 0.01);
        BOOST_CHECK_CLOSE(A.min()(1), -1.0, 0.01);
        BOOST_CHECK_CLOSE(A.min()(2), -1.0, 0.01);

        BOOST_CHECK_CLOSE(A.max()(0), 1.0, 0.01);
        BOOST_CHECK_CLOSE(A.max()(1), 1.0, 0.01);
        BOOST_CHECK_CLOSE(A.max()(2), 1.0, 0.01);
    }

    {
        EigenVector3<T> const min_coord = EigenVector3<T>(-1.0, -1.0, -1.0);
        EigenVector3<T> const max_coord = EigenVector3<T>(1.0, 1.0, 1.0);

        geometry::AABBEigen<T> A;

        A.max() = min_coord;
        A.min() = max_coord;

        BOOST_CHECK(geometry::is_valid(A) == false);

        BOOST_CHECK_CLOSE(A.max()(0), -1.0, 0.01);
        BOOST_CHECK_CLOSE(A.max()(1), -1.0, 0.01);
        BOOST_CHECK_CLOSE(A.max()(2), -1.0, 0.01);

        BOOST_CHECK_CLOSE(A.min()(0), 1.0, 0.01);
        BOOST_CHECK_CLOSE(A.min()(1), 1.0, 0.01);
        BOOST_CHECK_CLOSE(A.min()(2), 1.0, 0.01);
    }

    {
        EigenVector3<T> const min_coord = EigenVector3<T>(-1.0, -1.0, -1.0);
        EigenVector3<T> const max_coord = EigenVector3<T>(1.0, 1.0, 1.0);

        geometry::AABBEigen<T> A;

        A.min() = min_coord;
        A.max() = max_coord;

        geometry::AABBEigen<T> B = A;
        geometry::AABBEigen<T> C(A);

        BOOST_CHECK(geometry::is_valid(B) == true);
        BOOST_CHECK(geometry::is_valid(C) == true);

        BOOST_CHECK_CLOSE(B.min()(0), -1.0, 0.01);
        BOOST_CHECK_CLOSE(B.min()(1), -1.0, 0.01);
        BOOST_CHECK_CLOSE(B.min()(2), -1.0, 0.01);

        BOOST_CHECK_CLOSE(B.max()(0), 1.0, 0.01);
        BOOST_CHECK_CLOSE(B.max()(1), 1.0, 0.01);
        BOOST_CHECK_CLOSE(B.max()(2), 1.0, 0.01);

        BOOST_CHECK_CLOSE(C.min()(0), -1.0, 0.01);
        BOOST_CHECK_CLOSE(C.min()(1), -1.0, 0.01);
        BOOST_CHECK_CLOSE(C.min()(2), -1.0, 0.01);

        BOOST_CHECK_CLOSE(C.max()(0), 1.0, 0.01);
        BOOST_CHECK_CLOSE(C.max()(1), 1.0, 0.01);
        BOOST_CHECK_CLOSE(C.max()(2), 1.0, 0.01);
    }
}

BOOST_AUTO_TEST_SUITE_END();
