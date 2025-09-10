#include <geometry.h>
#include <eigenhelperall.h>

#define BOOST_AUTO_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include <boost/test/unit_test_suite.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/test/test_tools.hpp>

using T = double;
BOOST_AUTO_TEST_SUITE(geometry);

BOOST_AUTO_TEST_CASE(capsule_test)
{
    {
        geometry::CapsuleEigen<T> const A;

        BOOST_CHECK(geometry::is_valid(A) == true);

        BOOST_CHECK_CLOSE(A.point0()(0), 0.0, 0.01);
        BOOST_CHECK_CLOSE(A.point0()(1), 0.0, 0.01);
        BOOST_CHECK_CLOSE(A.point0()(2), 0.0, 0.01);

        BOOST_CHECK_CLOSE(A.point1()(0), 0.0, 0.01);
        BOOST_CHECK_CLOSE(A.point1()(1), 0.0, 0.01);
        BOOST_CHECK_CLOSE(A.point1()(2), 0.0, 0.01);

        BOOST_CHECK_CLOSE(A.radius(), 1.0, 0.01);
    }
    {
        EigenVector3<T> const point0 = EigenVector3<T>(1.0, 2.0, 3.0);
        EigenVector3<T> const point1 = EigenVector3<T>(4.0, 5.0, 6.0);
        T const radius = 7.0;

        geometry::CapsuleEigen<T> const A
            = geometry::make_capsule(radius, point0, point1);

        BOOST_CHECK(geometry::is_valid(A) == true);

        BOOST_CHECK_CLOSE(A.point0()(0), 1.0, 0.01);
        BOOST_CHECK_CLOSE(A.point0()(1), 2.0, 0.01);
        BOOST_CHECK_CLOSE(A.point0()(2), 3.0, 0.01);

        BOOST_CHECK_CLOSE(A.point1()(0), 4.0, 0.01);
        BOOST_CHECK_CLOSE(A.point1()(1), 5.0, 0.01);
        BOOST_CHECK_CLOSE(A.point1()(2), 6.0, 0.01);

        BOOST_CHECK_CLOSE(A.radius(), 7.0, 0.01);
    }
    {
        EigenVector3<T> const point0 = EigenVector3<T>(1.0, 2.0, 3.0);
        EigenVector3<T> const point1 = EigenVector3<T>(4.0, 5.0, 6.0);
        T const radius = 7.0;

        geometry::CapsuleEigen<T> A;

        A.point0() = point0;
        A.point1() = point1;
        A.set_radius(radius);

        BOOST_CHECK(geometry::is_valid(A) == true);

        BOOST_CHECK_CLOSE(A.point0()(0), 1.0, 0.01);
        BOOST_CHECK_CLOSE(A.point0()(1), 2.0, 0.01);
        BOOST_CHECK_CLOSE(A.point0()(2), 3.0, 0.01);

        BOOST_CHECK_CLOSE(A.point1()(0), 4.0, 0.01);
        BOOST_CHECK_CLOSE(A.point1()(1), 5.0, 0.01);
        BOOST_CHECK_CLOSE(A.point1()(2), 6.0, 0.01);

        BOOST_CHECK_CLOSE(A.radius(), 7.0, 0.01);
    }
    {
        EigenVector3<T> const point0 = EigenVector3<T>(1.0, 2.0, 3.0);
        EigenVector3<T> const point1 = EigenVector3<T>(4.0, 5.0, 6.0);
        T const radius = 7.0;

        geometry::CapsuleEigen<T> const A
            = geometry::make_capsule(radius, point0, point1);
        geometry::CapsuleEigen<T> const B(A);
        geometry::CapsuleEigen<T> const C = B;

        BOOST_CHECK(geometry::is_valid(B) == true);

        BOOST_CHECK_CLOSE(B.point0()(0), 1.0, 0.01);
        BOOST_CHECK_CLOSE(B.point0()(1), 2.0, 0.01);
        BOOST_CHECK_CLOSE(B.point0()(2), 3.0, 0.01);

        BOOST_CHECK_CLOSE(B.point1()(0), 4.0, 0.01);
        BOOST_CHECK_CLOSE(B.point1()(1), 5.0, 0.01);
        BOOST_CHECK_CLOSE(B.point1()(2), 6.0, 0.01);

        BOOST_CHECK_CLOSE(B.radius(), 7.0, 0.01);

        BOOST_CHECK_CLOSE(C.point0()(0), 1.0, 0.01);
        BOOST_CHECK_CLOSE(C.point0()(1), 2.0, 0.01);
        BOOST_CHECK_CLOSE(C.point0()(2), 3.0, 0.01);

        BOOST_CHECK_CLOSE(C.point1()(0), 4.0, 0.01);
        BOOST_CHECK_CLOSE(C.point1()(1), 5.0, 0.01);
        BOOST_CHECK_CLOSE(C.point1()(2), 6.0, 0.01);

        BOOST_CHECK_CLOSE(C.radius(), 7.0, 0.01);
    }

    {
        EigenVector3<T> const& center = EigenVector3<T>(0, 0, 0);
        EigenVector3<T> const& axis = EigenVector3<T>(0, 0, 1);
        T const& height = 2;
        const T radius = 3;

        geometry::CylinderEigen<T> const& C
            = geometry::make_cylinder(radius, height, axis, center);
        geometry::CapsuleEigen<T> const& A = geometry::convert(C);

        BOOST_CHECK_CLOSE(A.point0()(0), 0.0, 0.01);
        BOOST_CHECK_CLOSE(A.point0()(1), 0.0, 0.01);
        BOOST_CHECK_CLOSE(A.point0()(2), -1.0, 0.01);
        BOOST_CHECK_CLOSE(A.point1()(0), 0.0, 0.01);
        BOOST_CHECK_CLOSE(A.point1()(1), 0.0, 0.01);
        BOOST_CHECK_CLOSE(A.point1()(2), 1.0, 0.01);
        BOOST_CHECK_CLOSE(A.radius(), 3.0, 0.01);
    }
}

BOOST_AUTO_TEST_SUITE_END();
