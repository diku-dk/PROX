#include <geometry.h>
#include <tiny_math_types.h>

#define BOOST_AUTO_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include <boost/test/unit_test_suite.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/test/test_tools.hpp>

BOOST_AUTO_TEST_SUITE(geometry);

BOOST_AUTO_TEST_CASE(raycast_sphere)
{
    using std::sqrt;

    using T = double;

    EigenVector3<T> const center = EigenVector3<T>(0.0, 0.0, 0.0);
    T const radius = 1;

    geometry::Sphere<T> sphere = geometry::make_sphere(center, radius);

  // ray hitting straight on
    {
        EigenVector3<T> const r = EigenVector3<T>(0.0, 0.0, 1.0);
        EigenVector3<T> const p = EigenVector3<T>(0.0, 0.0, -3.0);

        geometry::RayEigen<T> const ray = geometry::make_ray(p, r);

        T length = 0;
        EigenVector3<T> q = EigenVector3<T>(0, 0, 0);

        bool hit = geometry::compute_raycast_sphere(ray, sphere, q, length);

        BOOST_CHECK(hit);

        BOOST_CHECK_CLOSE(length, 2.0, 0.01);
        BOOST_CHECK_CLOSE(q(0), 0.0, 0.01);
        BOOST_CHECK_CLOSE(q(1), 0.0, 0.01);
        BOOST_CHECK_CLOSE(q(2), -1.0, 0.01);
    }

  // ray shooting staright away
    {
        EigenVector3<T> const r = EigenVector3<T>(0.0, 0.0, -1.0);
        EigenVector3<T> const p = EigenVector3<T>(0.0, 0.0, -3.0);

        geometry::RayEigen<T> const ray = geometry::make_ray(p, r);

        T length = 0;
        EigenVector3<T> q = EigenVector3<T>(0, 0, 0);

        bool hit = geometry::compute_raycast_sphere(ray, sphere, q, length);

        BOOST_CHECK(!hit);
    }

  // ray starts inside
    {
        EigenVector3<T> const r = EigenVector3<T>(0.0, 0.0, -1.0);
        EigenVector3<T> const p = EigenVector3<T>(0.0, 0.0, 0.5);

        geometry::RayEigen<T> const ray = geometry::make_ray(p, r);

        T length = 0;
        EigenVector3<T> q = EigenVector3<T>(0, 0, 0);

        bool hit = geometry::compute_raycast_sphere(ray, sphere, q, length);

        BOOST_CHECK(!hit);
    }

  // ray shooting pass sphere
    {
        EigenVector3<T> const r = EigenVector3<T>(4.0, 0.0, 1.0);
        EigenVector3<T> const p = EigenVector3<T>(0.0, 0.0, -3.0);

        geometry::RayEigen<T> const ray = geometry::make_ray(p, r);

        T length = 0;
        EigenVector3<T> q = EigenVector3<T>(0, 0, 0);

        bool hit = geometry::compute_raycast_sphere(ray, sphere, q, length);

        BOOST_CHECK(!hit);
    }
}

BOOST_AUTO_TEST_SUITE_END();
