#include <geometry.h>
#include <tiny_math_types.h>

#define BOOST_AUTO_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include <boost/test/unit_test_suite.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/test/test_tools.hpp>

BOOST_AUTO_TEST_SUITE(geometry);

BOOST_AUTO_TEST_CASE(raycast_obb)
{
    using std::sqrt;

    using T = double;
    EigenVector3<T> const center = EigenVector3<T>(0.0, 0.0, 0.0);
    EigenVector3<T> const half_ext = EigenVector3<T>(1.0, 1.0, 1.0);
    EigenQuaternion<T> const q = EigenQuaternion<T>::Identity();

    geometry::OBBEigen<T> obb = geometry::make_obb<T>(center, q, half_ext);

  // Hit straigth on
    {
        EigenVector3<T> const r = EigenVector3<T>(0.0, 0.0, 1.0);
        EigenVector3<T> const p = EigenVector3<T>(0.0, 0.0, -3.0);

        geometry::RayEigen<T> const ray = geometry::make_ray(p, r);

        T length = 0;
        EigenVector3<T> q = EigenVector3<T>(0, 0, 0);

        bool hit = geometry::compute_raycast_obb(ray, obb, q, length);

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

        bool hit = geometry::compute_raycast_obb(ray, obb, q, length);

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

        bool hit = geometry::compute_raycast_obb(ray, obb, q, length);

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

        bool hit = geometry::compute_raycast_obb(ray, obb, q, length);

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

        bool hit = geometry::compute_raycast_obb(ray, obb, q, length);

        BOOST_CHECK(!hit);

        BOOST_CHECK_CLOSE(length, std::numeric_limits<T>::max(), 0.01);
        BOOST_CHECK_CLOSE(q(0), 0.0, 0.01);
        BOOST_CHECK_CLOSE(q(1), 0.0, 0.01);
        BOOST_CHECK_CLOSE(q(2), 0.0, 0.01);
    }

  // Obligue ray central hit
    {
        EigenVector3<T> const r = EigenVector3<T>(1.0, 1.0, 1.0);
        EigenVector3<T> const p = EigenVector3<T>(-1.0, -1.0, -2.0);

        geometry::RayEigen<T> const ray = geometry::make_ray(p, r);

        T length = 0;
        EigenVector3<T> q = EigenVector3<T>(0, 0, 0);

        bool hit = geometry::compute_raycast_obb(ray, obb, q, length);

        BOOST_CHECK(hit);

        BOOST_CHECK_CLOSE(length, 1.732050807568877, 0.01);
        BOOST_CHECK_CLOSE(q(0), 0.0, 0.01);
        BOOST_CHECK_CLOSE(q(1), 0.0, 0.01);
        BOOST_CHECK_CLOSE(q(2), -1.0, 0.01);
    }

  // Obligue ray corner hit
    {
        EigenVector3<T> const r = EigenVector3<T>(1.0, 1.0, 1.0);
        EigenVector3<T> const p = EigenVector3<T>(-1.0, -1.0, -2.9999);

        geometry::RayEigen<T> const ray = geometry::make_ray(p, r);

        T length = 0;
        EigenVector3<T> q = EigenVector3<T>(0, 0, 0);

        bool hit = geometry::compute_raycast_obb(ray, obb, q, length);

        BOOST_CHECK(hit);

        BOOST_CHECK_CLOSE(length, 2 * 1.732050807568877, 0.1);
        BOOST_CHECK_CLOSE(q(0), 1.0, 0.1);
        BOOST_CHECK_CLOSE(q(1), 1.0, 0.1);
        BOOST_CHECK_CLOSE(q(2), -1.0, 0.1);
    }

  // Obligue ray edge hit
    {
        EigenVector3<T> const r = EigenVector3<T>(1.0, 1.0, 1.0);
        EigenVector3<T> const p = EigenVector3<T>(0.0, -1.0, -2.0);

        geometry::RayEigen<T> const ray = geometry::make_ray(p, r);

        T length = 0;
        EigenVector3<T> q = EigenVector3<T>(0, 0, 0);

        bool hit = geometry::compute_raycast_obb(ray, obb, q, length);

        BOOST_CHECK(hit);

        BOOST_CHECK_CLOSE(length, 1.732050807568877, 0.01);
        BOOST_CHECK_CLOSE(q(0), 1.0, 0.01);
        BOOST_CHECK_CLOSE(q(1), 0.0, 0.01);
        BOOST_CHECK_CLOSE(q(2), -1.0, 0.01);
    }
}

BOOST_AUTO_TEST_CASE(raycast_strange_size_obb)
{
    using std::sqrt;

    using T = double;

    EigenVector3<T> const center = EigenVector3<T>(0.0, 0.0, 0.0);
    EigenVector3<T> const half_ext = EigenVector3<T>(1.0, 2.0, 3.0);
  //  Q const q        = Q::Rx(VT::pi_half() );
    EigenQuaternion<T> const q = EigenQuaternion<T>::Identity();

    geometry::OBBEigen<T> obb = geometry::make_obb<T>(center, q, half_ext);

    {
        EigenVector3<T> const r = EigenVector3<T>(2.0, 3.0, -1.0);
        EigenVector3<T> const p = EigenVector3<T>(-2.0, -3.0, 4.0);

        geometry::RayEigen<T> const ray = geometry::make_ray(p, r);

        T length = 0;
        EigenVector3<T> q = EigenVector3<T>(0, 0, 0);

        bool hit = geometry::compute_raycast_obb(ray, obb, q, length);

        BOOST_CHECK(hit);

        BOOST_CHECK_CLOSE(length, norm(r), 0.01);
        BOOST_CHECK_CLOSE(q(0), 0.0, 0.01);
        BOOST_CHECK_CLOSE(q(1), 0.0, 0.01);
        BOOST_CHECK_CLOSE(q(2), 3.0, 0.01);
    }
}

BOOST_AUTO_TEST_CASE(raycast_rotated_obb)
{
    using std::sqrt;

    using T = double;

    EigenVector3<T> const center = EigenVector3<T>(0.0, 0.0, 0.0);
    EigenVector3<T> const half_ext = EigenVector3<T>(1.0, 2.0, 3.0);
    EigenQuaternion<T> const q = Rotatex(std::numbers::pi_v<T> * 0.5f);

    geometry::OBBEigen<T> obb = geometry::make_obb<T>(center, q, half_ext);

    {
        EigenVector3<T> const r = EigenVector3<T>(2.0, 1.0, 3.0);
        EigenVector3<T> const p = EigenVector3<T>(-2.0, -4.0, -3.0);

        geometry::RayEigen<T> const ray = geometry::make_ray(p, r);

        T length = 0;
        EigenVector3<T> q = EigenVector3<T>(0, 0, 0);

        bool hit = geometry::compute_raycast_obb(ray, obb, q, length);

        BOOST_CHECK(hit);

        BOOST_CHECK_CLOSE(length, norm(r), 0.01);
        BOOST_CHECK_SMALL(q(0), 0.01);
        BOOST_CHECK_CLOSE(q(1), -3.0, 0.01);
        BOOST_CHECK_SMALL(q(2), 0.01);
    }
}

BOOST_AUTO_TEST_SUITE_END();
