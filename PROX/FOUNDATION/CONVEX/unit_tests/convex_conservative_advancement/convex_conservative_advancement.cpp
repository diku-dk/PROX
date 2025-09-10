#include <tiny.h>
#include <geometry.h>  // needed for geometry::Sphere
#include <convex.h>

#define BOOST_AUTO_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include <boost/test/unit_test_suite.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/test/test_tools.hpp>

#include <cmath>   // for std::fabs and std::sqrt

BOOST_AUTO_TEST_SUITE(conservative_advancement);

BOOST_AUTO_TEST_CASE(case_by_case_test)
{
    using std::fabs;

/*    typedef tiny::MathTypes<double> M;
    typedef M::value_traits VT;
    typedef M::vector3_type V;
    typedef M::real_type T;
    typedef M::quaternion_type Q;
    typedef M::coordsys_type C;*/
    using T = double;

    // Setup a central impact, two sphere hitting in perfect symmetry, totally independent of their rotationnal motion!
    {
        geometry::Sphere<T> const A;
        geometry::Sphere<T> const B;

        CoordSysEigen<T> T_a;
        CoordSysEigen<T> T_b;

        T_a.T() = EigenVector3<T>(-2.0, 0.0, 0.0);
        T_a.Q() = EigenQuaternion<T>::Identity();
        const EigenVector3<T> v_a = EigenVector3<T>(2.0, 0.0, 0.0);
        const EigenVector3<T> w_a = EigenVector3<T>(0.0, 5.0, 0.0);
        T const r_max_a = A.radius();

        T_b.T() = EigenVector3<T>(2.0, 0.0, 0.0);
        T_b.Q() = EigenQuaternion<T>::Identity();
        const EigenVector3<T> v_b = EigenVector3<T>(-2.0, 0.0, 0.0);
        const EigenVector3<T> w_b = EigenVector3<T>(0.0, 0.0, 5.0);
        T const r_max_b = B.radius();

        size_t const max_iterations = 100u;

        T const epsilon = 0.01;
        T toi = 0.0;
        size_t iterations;
        EigenVector3<T> p_a;
        EigenVector3<T> p_b;

        bool const impact = convex::conservative_advancement<T>(
            T_a, v_a, w_a, &A, r_max_a, T_b, v_b, w_b, &B, r_max_b, p_a, p_b,
            toi, iterations, epsilon, 1, max_iterations);

        BOOST_CHECK(impact);
        BOOST_CHECK_CLOSE(toi, 0.5, 1.0);

        BOOST_CHECK(fabs(p_a(0)) < epsilon);
        BOOST_CHECK(fabs(p_a(1)) < epsilon);
        BOOST_CHECK(fabs(p_a(2)) < epsilon);
        BOOST_CHECK(fabs(p_b(0)) < epsilon);
        BOOST_CHECK(fabs(p_b(1)) < epsilon);
        BOOST_CHECK(fabs(p_b(2)) < epsilon);
    }

    // Setup separating motion, two spheres moving away from each other
    {
        geometry::Sphere<T> const A;
        geometry::Sphere<T> const B;

        CoordSysEigen<T> T_a;
        CoordSysEigen<T> T_b;

        T_a.T() = EigenVector3<T>(-2.0, 0.0, 0.0);
        T_a.Q() = EigenQuaternion<T>::Identity();
        EigenVector3<T> const v_a = EigenVector3<T>(-2.0, 0.0, 0.0);
        EigenVector3<T> const w_a = EigenVector3<T>(0.0, 0.0, 0.0);
        T const r_max_a = A.radius();

        T_b.T() = EigenVector3<T>(2.0, 0.0, 0.0);
        T_b.Q() = EigenQuaternion<T>::Identity();
        EigenVector3<T> const v_b = EigenVector3<T>(2.0, 0.0, 0.0);
        EigenVector3<T> const w_b = EigenVector3<T>(0.0, 0.0, 0.0);
        T const r_max_b = B.radius();

        size_t const max_iterations = 100u;
        T const epsilon = 0.01;
        T toi = 0.0;
        size_t iterations;
        EigenVector3<T> p_a;
        EigenVector3<T> p_b;

        bool const impact = convex::conservative_advancement<T>(
            T_a, v_a, w_a, &A, r_max_a, T_b, v_b, w_b, &B, r_max_b, p_a, p_b,
            toi, iterations, epsilon, 1, max_iterations);

        BOOST_CHECK(!impact);
    }

    // Setup two sphere moving close by each other but never impacting
    {
        geometry::Sphere<T> const A;
        geometry::Sphere<T> const B;

        CoordSysEigen<T> T_a;
        CoordSysEigen<T> T_b;

        T_a.T() = EigenVector3<T>(-2.0, 1.01, 0.0);
        T_a.Q() = EigenQuaternion<T>::Identity();
        EigenVector3<T> const v_a = EigenVector3<T>(-2.0, 0.0, 0.0);
        EigenVector3<T> const w_a = EigenVector3<T>(0.0, 5.0, 0.0);
        T const r_max_a = A.radius();

        T_b.T() = EigenVector3<T>(2.0, -1.01, 0.0);
        T_b.Q() = EigenQuaternion<T>::Identity();
        EigenVector3<T> const v_b = EigenVector3<T>(2.0, 0.0, 0.0);
        EigenVector3<T> const w_b = EigenVector3<T>(0.0, 0.0, 5.0);
        T const r_max_b = B.radius();

        size_t const max_iterations = 100u;
        T const epsilon = 0.01;
        T toi = 0.0;
        size_t iterations;
        EigenVector3<T> p_a;
        EigenVector3<T> p_b;

        bool const impact = convex::conservative_advancement<T>(
            T_a, v_a, w_a, &A, r_max_a, T_b, v_b, w_b, &B, r_max_b, p_a, p_b,
            toi, iterations, epsilon, 1, max_iterations);

        BOOST_CHECK(!impact);
    }

    // Setup test case where only rotational motion is of importance!
    {
        using std::sqrt;

        convex::Cylinder<T> A;
        geometry::Sphere<T> B;

        A.half_height() = 10.0;
        A.radius() = 1.0;
        B.radius() = 1.0;

        CoordSysEigen<T> T_a;
        CoordSysEigen<T> T_b;

        T_a.T() = EigenVector3<T>(0.0, 0.0, 2.0);
        T_a.Q() = EigenQuaternion<T>::Identity();
        EigenVector3<T> const v_a = EigenVector3<T>(0.0, 0.0, 0.0);
        EigenVector3<T> const w_a
            = EigenVector3<T>(-std::numbers::pi_v<T>, 0.0, 0.0);
        T const r_max_a = sqrt(A.half_height() * A.half_height() + A.radius() * A.radius());

        T_b.T() = EigenVector3<T>(0.0, 10.0, 0.0);
        T_b.Q() = EigenQuaternion<T>::Identity();
        EigenVector3<T> const v_b = EigenVector3<T>(0.0, 0.0, 0.0);
        EigenVector3<T> const w_b = EigenVector3<T>(0.0, 0.0, 0.0);
        T const r_max_b = B.radius();

        size_t const max_iterations = 100u;
        T const epsilon = 0.01;
        T toi = 0.0;
        size_t iterations;
        EigenVector3<T> p_a;
        EigenVector3<T> p_b;

        bool const impact = convex::conservative_advancement<T>(
            T_a, v_a, w_a, &A, r_max_a, T_b, v_b, w_b, &B, r_max_b, p_a, p_b,
            toi, iterations, epsilon, 1, max_iterations);

        BOOST_CHECK(impact);
        BOOST_CHECK_CLOSE(toi, 0.5, 1.0);
        BOOST_CHECK(norm(p_a - p_b) < epsilon);
    }
}

BOOST_AUTO_TEST_SUITE_END();
