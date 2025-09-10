#include <tiny.h>
#include <geometry.h>  // needed for geometry::Box and geometry::Sphere
#include <convex.h>

#define BOOST_AUTO_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include <boost/test/unit_test_suite.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/test/test_tools.hpp>

#include <cmath>

/*using M = tiny::MathTypes<double>;
using Q = M::quaternion_type;
using V = M::vector3_type;
using T = M::real_type;
using C = M::coordsys_type;
using VT = M::value_traits;*/
using T = double;

BOOST_AUTO_TEST_SUITE(convex_compute_closest_points);

BOOST_AUTO_TEST_CASE(sphere_box_compile_test)
{
    {
        geometry::Sphere<T> B;
        geometry::BoxEigen<T> A;

        B.radius() = 0.5;
        A.half_extent() = EigenVector3<T>(50, .5, 50);

        CoordSysEigen<T> X_A;
        CoordSysEigen<T> X_B;

        EigenVector3<T> p_A = EigenVector3<T>(0, 0, 0);
        EigenVector3<T> p_B = EigenVector3<T>(0, 0, 0);

        X_A = CoordSysEigen<T>::identity();

        X_B.T() = EigenVector3<T>(0, 1.2, 0);
        X_B.Q() = EigenQuaternion<T>::Identity();

        convex::compute_closest_points<T>(X_A, &A, X_B, &B, p_A, p_B);
    }
    {

        geometry::BoxEigen<T> A;
        convex::Capsule<T> B;

        A.half_extent() = EigenVector3<T>(50, 1, 50);
        B.radius() = 0.5;

        CoordSysEigen<T> X_A;
        CoordSysEigen<T> X_B;

        EigenVector3<T> p_A = EigenVector3<T>(0, 0, 0);
        EigenVector3<T> p_B = EigenVector3<T>(0, 0, 0);

        X_A = CoordSysEigen<T>::identity();

        X_B.T() = EigenVector3<T>(0, 1.5, 0);
        X_B.Q() = EigenQuaternion<T>::Identity();

        convex::compute_closest_points<T>(X_A, &A, X_B, &B, p_A, p_B);
    }
}

BOOST_AUTO_TEST_CASE(sphere_sphere_testing)
{
    geometry::Sphere<T> const A;
    geometry::Sphere<T> const B;

    size_t const max_iterations = 100u;
    T const absolute_tolerance = T(10e-6);
    T const relative_tolerance = T(10e-6);
    T const stagnation_tolerance = T(10e-15);

  // Two unit-spheres placed ontop of each other
    {
        CoordSysEigen<T> X_A;
        CoordSysEigen<T> X_B;

        EigenVector3<T> p_A;
        EigenVector3<T> p_B;
        size_t iterations = 0u;
        size_t status = 0u;
        T distance = std::numeric_limits<T>::max();

        X_A = CoordSysEigen<T>::identity();
        X_B = CoordSysEigen<T>::identity();
        X_B.T() = EigenVector3<T>(0, 2, 0);

        convex::compute_closest_points<T>(
            X_A, &A, X_B, &B, p_A, p_B, distance, iterations, status,
            absolute_tolerance, relative_tolerance, stagnation_tolerance,
            max_iterations);

        T true_distance = norm(X_A.T() - X_B.T()) - 2;

        BOOST_CHECK(fabs(distance - true_distance) < 10e-6);

        BOOST_CHECK(p_A(0) == 0);
        BOOST_CHECK(p_A(1) == 1);
        BOOST_CHECK(p_A(2) == 0);

        BOOST_CHECK(p_B(0) == 0);
        BOOST_CHECK(p_B(1) == 1);
        BOOST_CHECK(p_B(2) == 0);
    }
  // Two unit-spheres overlapping but both placed on the x-axis
    {
        CoordSysEigen<T> X_A;
        CoordSysEigen<T> X_B;

        EigenVector3<T> p_A;
        EigenVector3<T> p_B;
        size_t iterations = 0u;
        size_t status = 0u;
        T distance = std::numeric_limits<T>::max();

        X_A = CoordSysEigen<T>::identity();
        X_B = CoordSysEigen<T>::identity();
        X_B.T()(0) = T(2.5);

        convex::compute_closest_points<T>(
            X_A, &A, X_B, &B, p_A, p_B, distance, iterations, status,
            absolute_tolerance, relative_tolerance, stagnation_tolerance,
            max_iterations);

        T true_distance = norm(X_A.T() - X_B.T()) - 2;

        BOOST_CHECK(fabs(distance - true_distance) < 10e-6);

        BOOST_CHECK(p_A(0) == 1);
        BOOST_CHECK(p_A(1) == 0);
        BOOST_CHECK(p_A(2) == 0);

        BOOST_CHECK(p_B(0) == T(1.5));
        BOOST_CHECK(p_B(1) == 0);
        BOOST_CHECK(p_B(2) == 0);
    }
  // Two unit-spheres exactly touching in one point (= one intersection point) and but both placed on the x-axis
    {
        CoordSysEigen<T> X_A;
        CoordSysEigen<T> X_B;

        EigenVector3<T> p_A;
        EigenVector3<T> p_B;
        size_t iterations = 0u;
        size_t status = 0u;
        T distance = std::numeric_limits<T>::max();

        X_A.T() = EigenVector3<T>(0, 0, 0);
        X_A.T()(0) = -2.0;
        X_A.Q() = EigenQuaternion<T>::Identity();
        X_B.T() = EigenVector3<T>(0, 0, 0);
        X_B.Q() = EigenQuaternion<T>::Identity();

        convex::compute_closest_points<T>(
            X_A, &A, X_B, &B, p_A, p_B, distance, iterations, status,
            absolute_tolerance, relative_tolerance, stagnation_tolerance,
            max_iterations);
    }
  // Two unit-spheres non-overlapping but both placed on the x-axis
    {
        CoordSysEigen<T> X_A;
        CoordSysEigen<T> X_B;

        EigenVector3<T> p_A;
        EigenVector3<T> p_B;
        size_t iterations = 0u;
        size_t status = 0u;
        T distance = std::numeric_limits<T>::max();

        X_A.T() = EigenVector3<T>(0, 0, 0);
        X_A.T()(0) = -2.5;
        X_A.Q() = EigenQuaternion<T>::Identity();
        X_B.T() = EigenVector3<T>(0, 0, 0);
        X_B.Q() = EigenQuaternion<T>::Identity();

        convex::compute_closest_points<T>(
            X_A, &A, X_B, &B, p_A, p_B, distance, iterations, status,
            absolute_tolerance, relative_tolerance, stagnation_tolerance,
            max_iterations);
    }
}

BOOST_AUTO_TEST_CASE(random_test)
{
    geometry::Sphere<T> const A;
    geometry::Sphere<T> const B;

    size_t const max_iterations = 100u;
    T const absolute_tolerance = T(10e-6);
    T const relative_tolerance = T(10e-10);
    T const stagnation_tolerance = T(0.0);

    for (size_t i = 0; i < 100u; ++i)
    {
        CoordSysEigen<T> X_A;
        CoordSysEigen<T> X_B;

        EigenVector3<T> p_A;
        EigenVector3<T> p_B;
        size_t iterations = 0u;
        size_t status = 0u;
        T distance = std::numeric_limits<T>::max();

        X_A.T() = randomEigen<T>(-2.0, 2.0);
        Eigen::Matrix<T, 4, 1> random_vals
            = Eigen::Matrix<T, 4, 1>::Random()
                  .cwiseAbs(); // Values between 0 and 1

        // Create quaternion from the random values (w, x, y, z)
        Eigen::Quaternion<T> q(random_vals(0), random_vals(1), random_vals(2),
                               random_vals(3));
        X_A.Q() = (q).normalized();
        X_B.T() = randomEigen<T>(-2.0, 2.0);
        random_vals = Eigen::Matrix<T, 4, 1>::Random()
                          .cwiseAbs(); // Values between 0 and 1

        // Create quaternion from the random values (w, x, y, z)
        q = Eigen::Quaternion<T>(random_vals(0), random_vals(1), random_vals(2),
                                 random_vals(3));
        X_B.Q() = (q).normalized();

        convex::compute_closest_points<T>(
            X_A, &A, X_B, &B, p_A, p_B, distance, iterations, status,
            absolute_tolerance, relative_tolerance, stagnation_tolerance,
            max_iterations);

        T const true_distance = norm(X_A.T() - X_B.T()) - 2;

        if (true_distance > absolute_tolerance)
        {
            BOOST_CHECK_CLOSE(true_distance, distance, 0.05);

            BOOST_CHECK(status != convex::ABSOLUTE_CONVERGENCE);            // Would indicate penetration
            BOOST_CHECK(status != convex::INTERSECTION);                    // Would indicate penetration
            BOOST_CHECK(status != convex::ITERATING);                       // Would indicate internal error in CONVEX
            BOOST_CHECK(status != convex::EXCEEDED_MAX_ITERATIONS_LIMIT);   // Would indicate internal error in CONVEX
            BOOST_CHECK(status != convex::NON_DESCEND_DIRECTION);           // Would indicate internal error in CONVEX

            //BOOST_CHECK( status != convex::RELATIVE_CONVERGENCE );          // Would indicate convergence to positive distance
            //BOOST_CHECK( status != convex::SIMPLEX_EXPANSION_FAILED );      // Would indicate convergence to positive distance
            //BOOST_CHECK( status != convex::STAGNATION );                    // Would indicate convergence to positive distance
            //BOOST_CHECK( status != convex::LOWER_ERROR_BOUND_CONVERGENCE ); // Would indicate convergence to positive distance
        }
        else
        {
            BOOST_CHECK(0.0 <= distance);
            BOOST_CHECK(distance <= absolute_tolerance);

            BOOST_CHECK(status != convex::STAGNATION); // Can only occur in case of positive distance
            BOOST_CHECK(status != convex::LOWER_ERROR_BOUND_CONVERGENCE); // Can only occur in case of positive distance
            BOOST_CHECK(status != convex::RELATIVE_CONVERGENCE); // Can only occur in case of positive distance
            BOOST_CHECK(status != convex::SIMPLEX_EXPANSION_FAILED); // Can only occur in case of positive distance
            BOOST_CHECK(status != convex::ITERATING); // Would indicate internal error in CONVEX
            BOOST_CHECK(status != convex::EXCEEDED_MAX_ITERATIONS_LIMIT); // Would indicate internal error in CONVEX
            BOOST_CHECK(status != convex::NON_DESCEND_DIRECTION); // Would indicate internal error in CONVEX

            //BOOST_CHECK( status != convex::ABSOLUTE_CONVERGENCE );         // Indicates penetration
            //BOOST_CHECK( status != convex::INTERSECTION );                 // Indicates penetration
        }
    }
}

BOOST_AUTO_TEST_SUITE_END();
