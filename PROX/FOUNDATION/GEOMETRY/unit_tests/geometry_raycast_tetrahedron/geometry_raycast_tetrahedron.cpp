#include <geometry.h>
#include <tiny_math_types.h>

#define BOOST_AUTO_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include <boost/test/unit_test_suite.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/test/test_tools.hpp>

BOOST_AUTO_TEST_SUITE(geometry);

BOOST_AUTO_TEST_CASE(raycast_tetrahedron)
{

    using std::sqrt;

    using T = double;

    EigenVector3<T> const p0 = EigenVector3<T>(0.0, 0.0, 0.0);
    EigenVector3<T> const p1 = EigenVector3<T>(1.0, 0.0, 0.0);
    EigenVector3<T> const p2 = EigenVector3<T>(0.0, 1.0, 0.0);
    EigenVector3<T> const p3 = EigenVector3<T>(0.0, 0.0, 1.0);

    geometry::TetrahedronEigen<T> const tetrahedron
        = geometry::make_tetrahedron(p0, p1, p2, p3);

  // Ray hitting on bottom of tetrahedron
    {
        EigenVector3<T> const p = EigenVector3<T>(0.0, 0.0, -1.0);
        EigenVector3<T> const r = EigenVector3<T>(0.2, 0.2, 1.0);
        geometry::RayEigen<T> const ray = geometry::make_ray(p, r);

        T length = 0;
        EigenVector3<T> q = EigenVector3<T>(0, 0, 0);

        bool const hit = geometry::compute_raycast_tetrahedron(ray, tetrahedron, q, length);

        BOOST_CHECK(hit);
        BOOST_CHECK_CLOSE(length, norm(r), 0.01);
        BOOST_CHECK_CLOSE(q(0), 0.2, 0.01);
        BOOST_CHECK_CLOSE(q(1), 0.2, 0.01);
        BOOST_CHECK_CLOSE(q(2), 0.0, 0.01);
    }
  // Ray missing tetrahedron
    {
        EigenVector3<T> const p = EigenVector3<T>(0.0, 0.0, -1.0);
        EigenVector3<T> const r = EigenVector3<T>(-0.2, -0.2, 1.0);
        geometry::RayEigen<T> const ray = geometry::make_ray(p, r);

        T length = 0;
        EigenVector3<T> q = EigenVector3<T>(0, 0, 0);

        bool const hit = geometry::compute_raycast_tetrahedron(ray, tetrahedron, q, length);

        BOOST_CHECK(!hit);
    }
  // Ray hitting on bottom of tetrahedron only bottom is surface
    {
        std::vector<bool> surface_map(4, false);

        surface_map[3] = true;

        EigenVector3<T> const p = EigenVector3<T>(0.0, 0.0, -1.0);
        EigenVector3<T> const r = EigenVector3<T>(0.2, 0.2, 1.0);
        geometry::RayEigen<T> const ray = geometry::make_ray(p, r);

        T length = 0;
        EigenVector3<T> q = EigenVector3<T>(0, 0, 0);

        bool const hit = geometry::compute_raycast_tetrahedron(ray, tetrahedron, q, length, surface_map);

        BOOST_CHECK(hit);
        BOOST_CHECK_CLOSE(length, norm(r), 0.01);
        BOOST_CHECK_CLOSE(q(0), 0.2, 0.01);
        BOOST_CHECK_CLOSE(q(1), 0.2, 0.01);
        BOOST_CHECK_CLOSE(q(2), 0.0, 0.01);
    }
  // Ray from above hitting bottom of tetrahedron only bottom is surface
    {
        std::vector<bool> surface_map(4, false);

        surface_map[3] = true;

        EigenVector3<T> const p = EigenVector3<T>(0.0, 0.0, 1.0);
        EigenVector3<T> const r = EigenVector3<T>(0.2, 0.2, -1.0);
        geometry::RayEigen<T> const ray = geometry::make_ray(p, r);

        T length = 0;
        EigenVector3<T> q = EigenVector3<T>(0, 0, 0);

        bool const hit = geometry::compute_raycast_tetrahedron(ray, tetrahedron, q, length, surface_map);

        BOOST_CHECK(!hit);
    }

  // Stress testing... generating a bunch of random rays that all are hitting the oblique top plane
    for (unsigned int samples = 0u; samples < 1000u; ++samples)
    {
        EigenVector3<T> const noise = randomEigen<T>(0.0, 1.0);
        T const v1 = noise(0);
        T const v2 = noise(1);
        T const v3 = 1.0 - v1 - v2;
        T const w1 = (v3 < 0.0) ? v1 / (1.0 - v3) : v1;
        T const w2 = (v3 < 0.0) ? v2 / (1.0 - v3) : v2;
        T const w3 = (v3 < 0.0) ? 0.0 : v3;

        if (w3 == 0.0) // We do not wish to test for exact edge cases --- they are sensitive to finite precision errors
            continue;
        if (w3 == 0.0) // We do not wish to test for exact edge cases --- they are sensitive to finite precision errors
            continue;
        if (w3 == 0.0) // We do not wish to test for exact edge cases --- they are sensitive to finite precision errors
            continue;

        BOOST_CHECK_CLOSE((w1 + w2 + w3), 1.0, 0.01);

        BOOST_CHECK(w1 >= 0.0);
        BOOST_CHECK(w2 >= 0.0);
        BOOST_CHECK(w3 >= 0.0);

        BOOST_CHECK(w1 <= 1.0);
        BOOST_CHECK(w2 <= 1.0);
        BOOST_CHECK(w3 <= 1.0);

        EigenVector3<T> const hit_point = w1 * p1 + w2 * p2 + w3 * p3;
        EigenVector3<T> const ray_origin = noise * 10.0;
        EigenVector3<T> const ray_direction = hit_point - ray_origin;
        geometry::RayEigen<T> const ray
            = geometry::make_ray(ray_origin, ray_direction);

        T length = 0;
        EigenVector3<T> q = EigenVector3<T>(0, 0, 0);

        bool const hit = geometry::compute_raycast_tetrahedron(ray, tetrahedron, q, length);

        BOOST_CHECK(hit);
        BOOST_CHECK_CLOSE(length, norm(ray_direction), 0.01);

        for (unsigned int i = 0; i < 3u; ++i)
        {
            if (fabs(hit_point(i)) < std::numeric_limits<T>::epsilon() * 10)
                BOOST_CHECK_SMALL(q(i), std::numeric_limits<T>::epsilon() * 10);
            else
                BOOST_CHECK_CLOSE(q(i), hit_point(i), 0.01);
        }
    }
}

BOOST_AUTO_TEST_SUITE_END();
