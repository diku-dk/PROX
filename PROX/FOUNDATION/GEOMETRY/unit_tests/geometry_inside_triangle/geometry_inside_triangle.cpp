#include <geometry.h>
#include <eigenhelperall.h>

#define BOOST_AUTO_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include <boost/test/unit_test_suite.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/test/test_tools.hpp>

using T = float;

BOOST_AUTO_TEST_SUITE(geometry);

BOOST_AUTO_TEST_CASE(inside_triangle_test)
{
    EigenVector3<T> const p0 = EigenVector3<T>(0.0, 0.0, 0.0);
    EigenVector3<T> const p1 = EigenVector3<T>(1.0, 0.0, 0.0);
    EigenVector3<T> const p2 = EigenVector3<T>(0.0, 1.0, 0.0);

    geometry::Triangle<T> const A = geometry::make_triangle(p0, p1, p2);

    BOOST_CHECK(geometry::is_valid(A) == true);

  // Close to corners are allways inside
    {
        EigenVector3<T> const k0 = EigenVector3<T>(0.0001, 0.0001, 0.0);
        EigenVector3<T> const k1 = EigenVector3<T>(0.9999, 0.0, 0.0);
        EigenVector3<T> const k2 = EigenVector3<T>(0.0, 0.9999, 0.0);

        bool test0 = geometry::inside_triangle(k0, A, true);
        bool test1 = geometry::inside_triangle(k1, A, true);
        bool test2 = geometry::inside_triangle(k2, A, true);

        BOOST_CHECK(test0);
        BOOST_CHECK(test1);
        BOOST_CHECK(test2);

        test0 = geometry::inside_triangle(k0, A, false);
        test1 = geometry::inside_triangle(k1, A, false);
        test2 = geometry::inside_triangle(k2, A, false);

        BOOST_CHECK(test0);
        BOOST_CHECK(test1);
        BOOST_CHECK(test2);
    }
  // Midpoint is allways inside
    {
        bool const test0
            = geometry::inside_triangle<T>((p0 + p1 + p2) / 3.0, A, true);
        bool const test1
            = geometry::inside_triangle<T>((p0 + p1 + p2) / 3.0, A, false);

        BOOST_CHECK(test0);
        BOOST_CHECK(test1);
    }
  // Some arbitary points
    {
        EigenVector3<T> const k0
            = EigenVector3<T>(0.2, 0.2, 0.0);  // inside and on plane
        EigenVector3<T> const k1
            = EigenVector3<T>(0.2, 0.2, 0.2); // inside and above plane
        EigenVector3<T> const k2
            = EigenVector3<T>(0.2, 0.2, -0.2); // inside and below plane
        EigenVector3<T> const k3
            = EigenVector3<T>(-1.0, -1.0, 0.0); // outside and on plane
        EigenVector3<T> const k4
            = EigenVector3<T>(-1.0, -1.0, 0.2); // outside and above plane
        EigenVector3<T> const k5
            = EigenVector3<T>(-1.0, -1.0, -0.2); // outside and below plane

        bool test0 = geometry::inside_triangle(k0, A, true);
        bool test1 = geometry::inside_triangle(k1, A, true);
        bool test2 = geometry::inside_triangle(k2, A, true);
        bool test3 = geometry::inside_triangle(k3, A, true);
        bool test4 = geometry::inside_triangle(k4, A, true);
        bool test5 = geometry::inside_triangle(k5, A, true);

        BOOST_CHECK(test0);
        BOOST_CHECK(!test1);
        BOOST_CHECK(!test2);
        BOOST_CHECK(!test3);
        BOOST_CHECK(!test4);
        BOOST_CHECK(!test5);

        test0 = geometry::inside_triangle(k0, A, false);
        test1 = geometry::inside_triangle(k1, A, false);
        test2 = geometry::inside_triangle(k2, A, false);
        test3 = geometry::inside_triangle(k3, A, false);
        test4 = geometry::inside_triangle(k4, A, false);
        test5 = geometry::inside_triangle(k5, A, false);

        BOOST_CHECK(test0);
        BOOST_CHECK(test1);
        BOOST_CHECK(test2);
        BOOST_CHECK(!test3);
        BOOST_CHECK(!test4);
        BOOST_CHECK(!test5);
    }
  // Different outside tests -- trying to generate all different cases....
    {
        EigenVector3<T> const k0 = EigenVector3<T>(0.5, -0.1, 0.0);
        EigenVector3<T> const k1 = EigenVector3<T>(1.0, 1.0, 0.0);
        EigenVector3<T> const k2 = EigenVector3<T>(-0.1, 0.5, 0.0);
        EigenVector3<T> const k3 = EigenVector3<T>(-1.0, -1.0, 0.0);
        EigenVector3<T> const k4 = EigenVector3<T>(0.1, 1.0, 0.0);
        EigenVector3<T> const k5 = EigenVector3<T>(1.0, 0.1, 0.0);

        bool const test0 = geometry::inside_triangle(k0, A, false);
        bool const test1 = geometry::inside_triangle(k1, A, false);
        bool const test2 = geometry::inside_triangle(k2, A, false);
        bool const test3 = geometry::inside_triangle(k3, A, false);
        bool const test4 = geometry::inside_triangle(k4, A, false);
        bool const test5 = geometry::inside_triangle(k5, A, false);

        BOOST_CHECK(!test0);
        BOOST_CHECK(!test1);
        BOOST_CHECK(!test2);
        BOOST_CHECK(!test3);
        BOOST_CHECK(!test4);
        BOOST_CHECK(!test5);
    }
}

BOOST_AUTO_TEST_SUITE_END();
