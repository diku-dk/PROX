#include <geometry.h>
#include <tiny_math_types.h>

#define BOOST_AUTO_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include <boost/test/unit_test_suite.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/test/test_tools.hpp>

BOOST_AUTO_TEST_SUITE(geometry);

BOOST_AUTO_TEST_CASE(raycast_plane)
{

    using std::sqrt;

    using T = double;

    EigenVector3<T> const normal = EigenVector3<T>(1.0, 1.0, 1.0);
    T const offset = norm(normal);
    geometry::Plane<T> const plane = geometry::make_plane(unit(normal), offset);

    BOOST_CHECK(geometry::is_valid(plane));

  // Orthogonal hit from back-side
    {
        EigenVector3<T> const p = EigenVector3<T>(0.0, 0.0, 0.0);
        EigenVector3<T> const r = EigenVector3<T>(1.0, 1.0, 1.0);
        geometry::RayEigen<T> const ray = geometry::make_ray(p, r);

        T length = 0;
        EigenVector3<T> q = EigenVector3<T>(0, 0, 0);

        bool hit = geometry::compute_raycast_plane(ray, plane, q, length);

        BOOST_CHECK(hit);

        BOOST_CHECK_CLOSE(length, offset, 0.01);
        BOOST_CHECK_CLOSE(q(0), 1.0, 0.01);
        BOOST_CHECK_CLOSE(q(1), 1.0, 0.01);
        BOOST_CHECK_CLOSE(q(2), 1.0, 0.01);
    }

  // Orthogonal hit from back-side with front-face only
    {
        EigenVector3<T> const p = EigenVector3<T>(0.0, 0.0, 0.0);
        EigenVector3<T> const r = EigenVector3<T>(1.0, 1.0, 1.0);
        geometry::RayEigen<T> const ray = geometry::make_ray(p, r);

        T length = 0;
        EigenVector3<T> q = EigenVector3<T>(0, 0, 0);

        bool hit = geometry::compute_raycast_plane(ray, plane, q, length, true);

        BOOST_CHECK(!hit);
    }

  // Orthogonal hit from front-side but with ray origin on back-side
    {
        EigenVector3<T> const p = EigenVector3<T>(0.0, 0.0, 0.0);
        EigenVector3<T> const r = EigenVector3<T>(-1.0, -1.0, -1.0);
        geometry::RayEigen<T> const ray = geometry::make_ray(p, r);

        T length = 0;
        EigenVector3<T> q = EigenVector3<T>(0, 0, 0);

        bool hit = geometry::compute_raycast_plane(ray, plane, q, length, true);

        BOOST_CHECK(!hit);
    }

  // Parallel ray and plane
    {
        EigenVector3<T> const p = EigenVector3<T>(0.0, 0.0, 0.0);
        EigenVector3<T> const r = EigenVector3<T>(-1.0, 1.0, 0.0);
        geometry::RayEigen<T> const ray = geometry::make_ray(p, r);

        T length = 0;
        EigenVector3<T> q = EigenVector3<T>(0, 0, 0);

        bool hit = geometry::compute_raycast_plane(ray, plane, q, length);

        BOOST_CHECK(!hit);
    }

  // Oblique hit
    {
        EigenVector3<T> const p = EigenVector3<T>(2.0, 1.0, 1.0);
        EigenVector3<T> const r = EigenVector3<T>(-1.0, 0.0, 0.0);
        geometry::RayEigen<T> const ray = geometry::make_ray(p, r);

        T length = 0;
        EigenVector3<T> q = EigenVector3<T>(0, 0, 0);

        bool hit = geometry::compute_raycast_plane(ray, plane, q, length);

        BOOST_CHECK(hit);

        BOOST_CHECK_CLOSE(length, 1.0, 0.01);
        BOOST_CHECK_CLOSE(q(0), 1.0, 0.01);
        BOOST_CHECK_CLOSE(q(1), 1.0, 0.01);
        BOOST_CHECK_CLOSE(q(2), 1.0, 0.01);
    }
}

BOOST_AUTO_TEST_SUITE_END();
