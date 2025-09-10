#include <geometry.h>
#include <eigenhelperall.h>

#define BOOST_AUTO_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include <boost/test/unit_test_suite.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/test/test_tools.hpp>

BOOST_AUTO_TEST_SUITE(geometry);

BOOST_AUTO_TEST_CASE(raycast_aabb)
{
    using std::sqrt;

    using T = float;

    EigenVector3<T> const min_coord = EigenVector3<T>(-1.0, -1.0, -1.0);
    EigenVector3<T> const max_coord = EigenVector3<T>(1.0, 1.0, 1.0);

    geometry::AABBEigen<T> aabb = geometry::make_aabb(min_coord, max_coord);

    BOOST_CHECK(geometry::is_valid(aabb));

  // Hit straigth on
    {
        EigenVector3<T> const r = EigenVector3<T>(0.0, 0.0, 1.0);
        EigenVector3<T> const p = EigenVector3<T>(0.0, 0.0, -3.0);

        geometry::RayEigen<T> const ray = geometry::make_ray(p, r);

        T length = 0;
        EigenVector3<T> q = EigenVector3<T>(0, 0, 0);

        bool hit = geometry::compute_raycast_aabb(ray, aabb, q, length);

        BOOST_CHECK(hit);

        BOOST_CHECK_CLOSE(length, 2.0, 0.01);
        BOOST_CHECK_CLOSE(q(0), 0.0, 0.01);
        BOOST_CHECK_CLOSE(q(1), 0.0, 0.01);
        BOOST_CHECK_CLOSE(q(2), -1.0, 0.01);
    }

  // Hit straigth on corner
    {
        EigenVector3<T> const r = EigenVector3<T>(0.0, 0.0, 1.0);
        EigenVector3<T> const p = EigenVector3<T>(-1.0, -1.0, -3.0);

        geometry::RayEigen<T> const ray = geometry::make_ray(p, r);

        T length = 0;
        EigenVector3<T> q = EigenVector3<T>(0, 0, 0);

        bool hit = geometry::compute_raycast_aabb(ray, aabb, q, length);

        BOOST_CHECK(hit);

        BOOST_CHECK_CLOSE(length, 2.0, 0.01);
        BOOST_CHECK_CLOSE(q(0), -1.0, 0.01);
        BOOST_CHECK_CLOSE(q(1), -1.0, 0.01);
        BOOST_CHECK_CLOSE(q(2), -1.0, 0.01);
    }

  // Hit straigth on edge
    {
        EigenVector3<T> const r = EigenVector3<T>(0.0, 0.0, 1.0);
        EigenVector3<T> const p = EigenVector3<T>(0.0, -1.0, -3.0);

        geometry::RayEigen<T> const ray = geometry::make_ray(p, r);

        T length = 0;
        EigenVector3<T> q = EigenVector3<T>(0, 0, 0);

        bool hit = geometry::compute_raycast_aabb(ray, aabb, q, length);

        BOOST_CHECK(hit);

        BOOST_CHECK_CLOSE(length, 2.0, 0.01);
        BOOST_CHECK_CLOSE(q(0), 0.0, 0.01);
        BOOST_CHECK_CLOSE(q(1), -1.0, 0.01);
        BOOST_CHECK_CLOSE(q(2), -1.0, 0.01);
    }

  // Aligned ray no hitting
    {
        EigenVector3<T> const r = EigenVector3<T>(0.0, 0.0, 1.0);
        EigenVector3<T> const p = EigenVector3<T>(-2.0, 0.0, -3.0);

        geometry::RayEigen<T> const ray = geometry::make_ray(p, r);

        T length = 0;
        EigenVector3<T> q = EigenVector3<T>(0, 0, 0);

        bool hit = geometry::compute_raycast_aabb(ray, aabb, q, length);

        BOOST_CHECK(!hit);

        BOOST_CHECK_CLOSE(length, std::numeric_limits<T>::max(), 0.01);
        BOOST_CHECK_CLOSE(q(0), 0.0, 0.01);
        BOOST_CHECK_CLOSE(q(1), 0.0, 0.01);
        BOOST_CHECK_CLOSE(q(2), 0.0, 0.01);
    }

  // Obligue ray no hitting
    {
        EigenVector3<T> const r = EigenVector3<T>(1.0, 1.0, 1.0);
        EigenVector3<T> const p = EigenVector3<T>(0.0, 0.0, -10.0);

        geometry::RayEigen<T> const ray = geometry::make_ray(p, r);

        T length = 0;
        EigenVector3<T> q = EigenVector3<T>(0, 0, 0);

        bool hit = geometry::compute_raycast_aabb(ray, aabb, q, length);

        BOOST_CHECK(!hit);

        BOOST_CHECK_CLOSE(length, std::numeric_limits<T>::max(), 0.01);
        BOOST_CHECK_CLOSE(q(0), 0.0, 0.01);
        BOOST_CHECK_CLOSE(q(1), 0.0, 0.01);
        BOOST_CHECK_CLOSE(q(2), 0.0, 0.01);
    }

  // Obligue ray central hit
    {
        EigenVector3<T> const r = EigenVector3<T>(1.0, 1.0, 1.0);
        EigenVector3<T> const p = EigenVector3<T>(-2.0, -2.0, -3.0);

        geometry::RayEigen<T> const ray = geometry::make_ray(p, r);

        T length = 0;
        EigenVector3<T> q = EigenVector3<T>(0, 0, 0);

        bool hit = geometry::compute_raycast_aabb(ray, aabb, q, length);

        BOOST_CHECK(hit);

        BOOST_CHECK_CLOSE(length, 2 * 1.732050807568877, 0.01);
        BOOST_CHECK_CLOSE(q(0), 0.0, 0.01);
        BOOST_CHECK_CLOSE(q(1), 0.0, 0.01);
        BOOST_CHECK_CLOSE(q(2), -1.0, 0.01);
    }

  // Obligue ray corner hit
    {
        EigenVector3<T> const r = EigenVector3<T>(1.0, 1.0, 1.0);
        EigenVector3<T> const p = EigenVector3<T>(0.0, 0.0, -2.0);

        geometry::RayEigen<T> const ray = geometry::make_ray(p, r);

        T length = 0;
        EigenVector3<T> q = EigenVector3<T>(0, 0, 0);

        bool hit = geometry::compute_raycast_aabb(ray, aabb, q, length);

        BOOST_CHECK(hit);

        BOOST_CHECK_CLOSE(length, 1.732050807568877, 0.01);
        BOOST_CHECK_CLOSE(q(0), 1.0, 0.01);
        BOOST_CHECK_CLOSE(q(1), 1.0, 0.01);
        BOOST_CHECK_CLOSE(q(2), -1.0, 0.01);
    }

  // Obligue ray edge hit
    {
        EigenVector3<T> const r = EigenVector3<T>(1.0, 1.0, 1.0);
        EigenVector3<T> const p = EigenVector3<T>(0.0, -1.0, -2.0);

        geometry::RayEigen<T> const ray = geometry::make_ray(p, r);

        T length = 0;
        EigenVector3<T> q = EigenVector3<T>(0, 0, 0);

        bool hit = geometry::compute_raycast_aabb(ray, aabb, q, length);

        BOOST_CHECK(hit);

        BOOST_CHECK_CLOSE(length, 1.732050807568877, 0.01);
        BOOST_CHECK_CLOSE(q(0), 1.0, 0.01);
        BOOST_CHECK_CLOSE(q(1), 0.0, 0.01);
        BOOST_CHECK_CLOSE(q(2), -1.0, 0.01);
    }
}

BOOST_AUTO_TEST_SUITE_END();
