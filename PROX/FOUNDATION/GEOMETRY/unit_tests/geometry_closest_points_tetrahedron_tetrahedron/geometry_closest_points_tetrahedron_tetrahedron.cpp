#include <geometry.h>
#include <tiny.h>

#define BOOST_AUTO_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include <boost/test/unit_test_suite.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/test/test_tools.hpp>

#include <vector>

using T = float;

BOOST_AUTO_TEST_SUITE(geometry);

BOOST_AUTO_TEST_CASE(closest_points_tetrahedron_tetrahedron_test)
{
  // Get V-V contact
    {
        EigenVector3<T> const a0 = EigenVector3<T>(0.0, 0.0, 0.0);
        EigenVector3<T> const a1 = EigenVector3<T>(1.0, 0.0, 0.0);
        EigenVector3<T> const a2 = EigenVector3<T>(0.0, 1.0, 0.0);
        EigenVector3<T> const a3 = EigenVector3<T>(0.0, 0.0, 1.0);

        EigenVector3<T> const offset = EigenVector3<T>(2.0, 0.0, 0.0);

        EigenVector3<T> const b0 = EigenVector3<T>(0.0, 0.0, 0.0) + offset;
        EigenVector3<T> const b1 = EigenVector3<T>(1.0, 0.0, 0.0) + offset;
        EigenVector3<T> const b2 = EigenVector3<T>(0.0, 1.0, 0.0) + offset;
        EigenVector3<T> const b3 = EigenVector3<T>(0.0, 0.0, 1.0) + offset;

        geometry::TetrahedronEigen<T> const A
            = geometry::make_tetrahedron(a0, a1, a2, a3);
        geometry::TetrahedronEigen<T> const B
            = geometry::make_tetrahedron(b0, b1, b2, b3);

        EigenVector3<T> a;
        EigenVector3<T> b;
        T distance;

        geometry::closest_points_tetrahedron_tetrahedron(A, B, a, b, distance);

        BOOST_CHECK_CLOSE(distance, 1.0, 0.01);
        BOOST_CHECK_CLOSE(a(0), 1.0, 0.01);
        BOOST_CHECK_CLOSE(a(1), 0.0, 0.01);
        BOOST_CHECK_CLOSE(a(2), 0.0, 0.01);
        BOOST_CHECK_CLOSE(b(0), 2.0, 0.01);
        BOOST_CHECK_CLOSE(b(1), 0.0, 0.01);
        BOOST_CHECK_CLOSE(b(2), 0.0, 0.01);
    }

  // Get V-F contact
    {
        EigenVector3<T> const a0 = EigenVector3<T>(0.0, 0.0, 0.0);
        EigenVector3<T> const a1 = EigenVector3<T>(1.0, 0.0, 0.0);
        EigenVector3<T> const a2 = EigenVector3<T>(0.0, 1.0, 0.0);
        EigenVector3<T> const a3 = EigenVector3<T>(0.0, 0.0, 1.0);

        EigenVector3<T> const offset = EigenVector3<T>(2.0, 2.0, 2.0);

        EigenVector3<T> const b0 = EigenVector3<T>(0.0, 0.0, 0.0) + offset;
        EigenVector3<T> const b1 = EigenVector3<T>(1.0, 0.0, 0.0) + offset;
        EigenVector3<T> const b2 = EigenVector3<T>(0.0, 1.0, 0.0) + offset;
        EigenVector3<T> const b3 = EigenVector3<T>(0.0, 0.0, 1.0) + offset;

        geometry::TetrahedronEigen<T> const A
            = geometry::make_tetrahedron(a0, a1, a2, a3);
        geometry::TetrahedronEigen<T> const B
            = geometry::make_tetrahedron(b0, b1, b2, b3);

        EigenVector3<T> a;
        EigenVector3<T> b;
        T distance;

        geometry::closest_points_tetrahedron_tetrahedron(A, B, a, b, distance);

        BOOST_CHECK_CLOSE(distance, 2.88675141, 0.01);
        BOOST_CHECK_CLOSE(a(0), 0.33333, 0.01);
        BOOST_CHECK_CLOSE(a(1), 0.33333, 0.01);
        BOOST_CHECK_CLOSE(a(2), 0.33333, 0.01);
        BOOST_CHECK_CLOSE(b(0), 2.0, 0.01);
        BOOST_CHECK_CLOSE(b(1), 2.0, 0.01);
        BOOST_CHECK_CLOSE(b(2), 2.0, 0.01);
    }

  // Get V-E contact
    {
        EigenVector3<T> const a0 = EigenVector3<T>(0.0, 0.0, 0.0);
        EigenVector3<T> const a1 = EigenVector3<T>(1.0, 0.0, 0.0);
        EigenVector3<T> const a2 = EigenVector3<T>(0.0, 1.0, 0.0);
        EigenVector3<T> const a3 = EigenVector3<T>(0.0, 0.0, 1.0);

        EigenVector3<T> const offset = EigenVector3<T>(0.5, 0.0, -2.0);

        EigenVector3<T> const b0 = EigenVector3<T>(0.0, 0.0, 0.0) + offset;
        EigenVector3<T> const b1 = EigenVector3<T>(1.0, 0.0, 0.0) + offset;
        EigenVector3<T> const b2 = EigenVector3<T>(0.0, 1.0, 0.0) + offset;
        EigenVector3<T> const b3 = EigenVector3<T>(0.0, 0.0, 1.0) + offset;

        geometry::TetrahedronEigen<T> const A
            = geometry::make_tetrahedron(a0, a1, a2, a3);
        geometry::TetrahedronEigen<T> const B
            = geometry::make_tetrahedron(b0, b1, b2, b3);

        EigenVector3<T> a;
        EigenVector3<T> b;
        T distance;

        geometry::closest_points_tetrahedron_tetrahedron(A, B, a, b, distance);

        BOOST_CHECK_CLOSE(distance, 1.0, 0.01);

        BOOST_CHECK_CLOSE(a(0), 0.5, 0.01);
        BOOST_CHECK_CLOSE(a(1), 0.0, 0.01);
        BOOST_CHECK_CLOSE(a(2), 0.0, 0.01);

        BOOST_CHECK_CLOSE(b(0), 0.5, 0.01);
        BOOST_CHECK_CLOSE(b(1), 0.0, 0.01);
        BOOST_CHECK_CLOSE(b(2), -1.0, 0.01);
    }

  // Get E-E contact
    {
        EigenVector3<T> const a0 = EigenVector3<T>(0.0, 0.0, 0.0);
        EigenVector3<T> const a1 = EigenVector3<T>(1.0, 0.0, 0.0);
        EigenVector3<T> const a2 = EigenVector3<T>(0.0, 1.0, 0.0);
        EigenVector3<T> const a3 = EigenVector3<T>(0.0, 0.0, 1.0);

        EigenVector3<T> const offset = EigenVector3<T>(2.0, 2.0, -0.5);

        EigenVector3<T> const b0 = EigenVector3<T>(0.0, 0.0, 0.0) + offset;
        EigenVector3<T> const b1 = EigenVector3<T>(1.0, 0.0, 0.0) + offset;
        EigenVector3<T> const b2 = EigenVector3<T>(0.0, 1.0, 0.0) + offset;
        EigenVector3<T> const b3 = EigenVector3<T>(0.0, 0.0, 1.0) + offset;

        geometry::TetrahedronEigen<T> const A
            = geometry::make_tetrahedron(a0, a1, a2, a3);
        geometry::TetrahedronEigen<T> const B
            = geometry::make_tetrahedron(b0, b1, b2, b3);

        EigenVector3<T> a;
        EigenVector3<T> b;
        T distance;

        geometry::closest_points_tetrahedron_tetrahedron(A, B, a, b, distance);

        BOOST_CHECK_CLOSE(distance, 2.121320343559642, 0.01);

        BOOST_CHECK_CLOSE(a(0), 0.5, 0.01);
        BOOST_CHECK_CLOSE(a(1), 0.5, 0.01);
        BOOST_CHECK_CLOSE(a(2), 0.0, 0.01);

        BOOST_CHECK_CLOSE(b(0), 2.0, 0.01);
        BOOST_CHECK_CLOSE(b(1), 2.0, 0.01);
        BOOST_CHECK_CLOSE(b(2), 0.0, 0.01);
    }
}

BOOST_AUTO_TEST_SUITE_END();
