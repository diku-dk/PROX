#include <geometry.h>
#include <tiny_math_types.h>

#define BOOST_AUTO_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include <boost/test/unit_test_suite.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/test/test_tools.hpp>

BOOST_AUTO_TEST_SUITE(geometry);

BOOST_AUTO_TEST_CASE(raycast_triangle)
{

    using std::sqrt;

    using T = double;

    EigenVector3<T> const p0 = EigenVector3<T>(0.0, 0.0, 0.0);
    EigenVector3<T> const p1 = EigenVector3<T>(1.0, 0.0, 0.0);
    EigenVector3<T> const p2 = EigenVector3<T>(0.0, 1.0, 0.0);

    geometry::Triangle<T> const triangle = geometry::make_triangle(p0, p1, p2);

  // Ray hitting from front
    {
        EigenVector3<T> const p = EigenVector3<T>(0.0, 0.0, 1.0);
        EigenVector3<T> const r = EigenVector3<T>(0.2, 0.2, -1.0);
        geometry::RayEigen<T> const ray = geometry::make_ray(p, r);

        T length = 0;
        EigenVector3<T> q = EigenVector3<T>(0, 0, 0);

        bool hit = geometry::compute_raycast_triangle(ray, triangle, q, length, true);

        BOOST_CHECK(hit);
        BOOST_CHECK_CLOSE(length, norm(r), 0.01);
        BOOST_CHECK_CLOSE(q(0), 0.2, 0.01);
        BOOST_CHECK_CLOSE(q(1), 0.2, 0.01);
        BOOST_CHECK_CLOSE(q(2), 0.0, 0.01);
    }
  // Ray hitting from front
    {
        EigenVector3<T> const p = EigenVector3<T>(0.0, 0.0, 1.0);
        EigenVector3<T> const r = EigenVector3<T>(0.2, 0.2, -1.0);
        geometry::RayEigen<T> const ray = geometry::make_ray(p, r);

        T length = 0;
        EigenVector3<T> q = EigenVector3<T>(0, 0, 0);

        bool hit = geometry::compute_raycast_triangle(ray, triangle, q, length, false);

        BOOST_CHECK(hit);
        BOOST_CHECK_CLOSE(length, norm(r), 0.01);
        BOOST_CHECK_CLOSE(q(0), 0.2, 0.01);
        BOOST_CHECK_CLOSE(q(1), 0.2, 0.01);
        BOOST_CHECK_CLOSE(q(2), 0.0, 0.01);
    }
  // Ray hitting from back
    {
        EigenVector3<T> const p = EigenVector3<T>(0.0, 0.0, -1.0);
        EigenVector3<T> const r = EigenVector3<T>(0.2, 0.2, 1.0);
        geometry::RayEigen<T> const ray = geometry::make_ray(p, r);

        T length = 0;
        EigenVector3<T> q = EigenVector3<T>(0, 0, 0);

        bool hit = geometry::compute_raycast_triangle(ray, triangle, q, length, true);

        BOOST_CHECK(!hit);
    }
  // Ray hitting from back
    {
        EigenVector3<T> const p = EigenVector3<T>(0.0, 0.0, -1.0);
        EigenVector3<T> const r = EigenVector3<T>(0.2, 0.2, 1.0);
        geometry::RayEigen<T> const ray = geometry::make_ray(p, r);

        T length = 0;
        EigenVector3<T> q = EigenVector3<T>(0, 0, 0);

        bool hit = geometry::compute_raycast_triangle(ray, triangle, q, length, false);

        BOOST_CHECK(hit);
        BOOST_CHECK_CLOSE(length, norm(r), 0.01);
        BOOST_CHECK_CLOSE(q(0), 0.2, 0.01);
        BOOST_CHECK_CLOSE(q(1), 0.2, 0.01);
        BOOST_CHECK_CLOSE(q(2), 0.0, 0.01);
    }

  // Rays missing
    {
        EigenVector3<T> const p = EigenVector3<T>(0.0, 0.0, 1.0);
        EigenVector3<T> const r = EigenVector3<T>(0.6, 0.6, -1.0);
        geometry::RayEigen<T> const ray = geometry::make_ray(p, r);

        T length = 0;
        EigenVector3<T> q = EigenVector3<T>(0, 0, 0);

        bool hit = geometry::compute_raycast_triangle(ray, triangle, q, length, false);

        BOOST_CHECK(!hit);
    }
    {
        EigenVector3<T> const p = EigenVector3<T>(0.0, 0.0, 1.0);
        EigenVector3<T> const r = EigenVector3<T>(0.5, -0.1, -1.0);
        geometry::RayEigen<T> const ray = geometry::make_ray(p, r);

        T length = 0;
        EigenVector3<T> q = EigenVector3<T>(0, 0, 0);

        bool hit = geometry::compute_raycast_triangle(ray, triangle, q, length, false);

        BOOST_CHECK(!hit);
    }

  // Ray gracing
    {
        EigenVector3<T> const p = EigenVector3<T>(0.0, 0.0, -1.0);
        EigenVector3<T> const r = EigenVector3<T>(0.4999, 0.4999, 1.0);
        geometry::RayEigen<T> const ray = geometry::make_ray(p, r);

        T length = 0;
        EigenVector3<T> q = EigenVector3<T>(0, 0, 0);

        bool hit = geometry::compute_raycast_triangle(ray, triangle, q, length, false);

        BOOST_CHECK(hit);
        BOOST_CHECK_CLOSE(length, norm(r), 0.01);
        BOOST_CHECK_CLOSE(q(0), r(0), 0.01);
        BOOST_CHECK_CLOSE(q(1), r(1), 0.01);
        BOOST_CHECK_CLOSE(q(2), 0.0, 0.01);
    }
}

BOOST_AUTO_TEST_SUITE_END();
