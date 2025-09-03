#include <geometry.h>
#include <tiny.h>

#define BOOST_AUTO_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include <boost/test/unit_test_suite.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/test/test_tools.hpp>

using MT = tiny::MathTypes<float>;
using V = MT::vector3_type;
using T = MT::real_type;

BOOST_AUTO_TEST_SUITE(geometry);

BOOST_AUTO_TEST_CASE(sphere_test)
{

    {
        geometry::Sphere<V> A;
        BOOST_CHECK(geometry::is_valid(A) == true);

        BOOST_CHECK_CLOSE(A.center()(0), 0.0, 0.01);
        BOOST_CHECK_CLOSE(A.center()(1), 0.0, 0.01);
        BOOST_CHECK_CLOSE(A.center()(2), 0.0, 0.01);
        BOOST_CHECK_CLOSE(A.radius(), 1.0, 0.01);
    }

    {
        V const center = V::make(1.0, 2.0, 3.0);
        T const radius = 4.0;

        geometry::Sphere<V> A = geometry::make_sphere(center, radius);
        BOOST_CHECK(geometry::is_valid(A) == true);

        BOOST_CHECK_CLOSE(A.center()(0), 1.0, 0.01);
        BOOST_CHECK_CLOSE(A.center()(1), 2.0, 0.01);
        BOOST_CHECK_CLOSE(A.center()(2), 3.0, 0.01);
        BOOST_CHECK_CLOSE(A.radius(), 4.0, 0.01);
    }

    {
        V const center = V::make(1.0, 2.0, 3.0);
        T const radius = 4.0;

        geometry::Sphere<V> A;

        A.center() = center;
        A.radius() = radius;

        BOOST_CHECK(geometry::is_valid(A) == true);

        BOOST_CHECK_CLOSE(A.center()(0), 1.0, 0.01);
        BOOST_CHECK_CLOSE(A.center()(1), 2.0, 0.01);
        BOOST_CHECK_CLOSE(A.center()(2), 3.0, 0.01);
        BOOST_CHECK_CLOSE(A.radius(), 4.0, 0.01);
    }

    {
        V const center = V::make(1.0, 2.0, 3.0);
        T const radius = 4.0;

        geometry::Sphere<V> A = geometry::make_sphere(center, radius);

        geometry::Sphere<V> B = A;
        geometry::Sphere<V> C(A);

        BOOST_CHECK(geometry::is_valid(B) == true);
        BOOST_CHECK(geometry::is_valid(C) == true);

        BOOST_CHECK_CLOSE(B.center()(0), 1.0, 0.01);
        BOOST_CHECK_CLOSE(B.center()(1), 2.0, 0.01);
        BOOST_CHECK_CLOSE(B.center()(2), 3.0, 0.01);
        BOOST_CHECK_CLOSE(B.radius(), 4.0, 0.01);

        BOOST_CHECK_CLOSE(C.center()(0), 1.0, 0.01);
        BOOST_CHECK_CLOSE(C.center()(1), 2.0, 0.01);
        BOOST_CHECK_CLOSE(C.center()(2), 3.0, 0.01);
        BOOST_CHECK_CLOSE(C.radius(), 4.0, 0.01);
    }
}

BOOST_AUTO_TEST_SUITE_END();
