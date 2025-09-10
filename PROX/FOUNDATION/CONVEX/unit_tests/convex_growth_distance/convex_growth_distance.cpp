#include <tiny.h>
#include <geometry.h>  // needed for geometry shape types
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

BOOST_AUTO_TEST_SUITE(convex_growth_distance);

BOOST_AUTO_TEST_CASE(growth_distance_spheres)
{
    geometry::Sphere<T> const A;
    geometry::Sphere<T> const B;

    CoordSysEigen<T> X_A, X_B;
    EigenVector3<T> p_A, p_B;
    T growth_scale;
    size_t iterations;

    size_t const max_its = 100u;
    T const epsilon = 0.001;

  //--- separated spheres -------------------------------------------------------
    X_A.T() = EigenVector3<T>(1.0, 1.0, 0.0);
    X_B.T() = EigenVector3<T>(4.0, 1.0, 0.0);

    bool succes = convex::growth_distance<T>(
        X_A, &A, X_B, &B, p_A, p_B, growth_scale, iterations, epsilon, max_its);

    BOOST_CHECK(succes);
    BOOST_CHECK_CLOSE(growth_scale, 1.50, 1.0);
    BOOST_CHECK_CLOSE(p_A(0), 2.5, 1.0);
    BOOST_CHECK_CLOSE(p_A(1), 1.0, 0.01);
    BOOST_CHECK(p_A(2) == 0);
    BOOST_CHECK_CLOSE(p_B(0), 2.5, 1.0);
    BOOST_CHECK_CLOSE(p_B(1), 1.0, 0.01);
    BOOST_CHECK(p_B(2) == 0);

  //--- touching spheres -------------------------------------------------------
    X_A.T() = EigenVector3<T>(1.0, 1.0, 0.0);
    X_B.T() = EigenVector3<T>(3.0, 1.0, 0.0);

    succes = convex::growth_distance<T>(
        X_A, &A, X_B, &B, p_A, p_B, growth_scale, iterations, epsilon, max_its);

    BOOST_CHECK(succes);
    BOOST_CHECK_CLOSE(growth_scale, 1.00, 1.0);
    BOOST_CHECK_CLOSE(p_A(0), 2.0, 1.0);
    BOOST_CHECK_CLOSE(p_A(1), 1.0, 0.01);
    BOOST_CHECK(p_A(2) == 0);
    BOOST_CHECK_CLOSE(p_B(0), 2.0, 1.0);
    BOOST_CHECK_CLOSE(p_B(1), 1.0, 0.01);
    BOOST_CHECK(p_B(2) == 0);

  //---- slight overlapping spheres a long x-axis ------------------------------
    X_A.T() = EigenVector3<T>(1.0, 1.0, 0.0);
    X_B.T() = EigenVector3<T>(2.5, 1.0, 0.0);

    succes = convex::growth_distance<T>(
        X_A, &A, X_B, &B, p_A, p_B, growth_scale, iterations, epsilon, max_its);

    BOOST_CHECK(succes);
    BOOST_CHECK_CLOSE(growth_scale, 0.75, 1.0);
    BOOST_CHECK_CLOSE(p_A(0), 1.75, 1.0);
    BOOST_CHECK_CLOSE(p_A(1), 1.0, 0.01);
    BOOST_CHECK(p_A(2) == 0);
    BOOST_CHECK_CLOSE(p_B(0), 1.75, 1.0);
    BOOST_CHECK_CLOSE(p_B(1), 1.0, 0.01);
    BOOST_CHECK(p_B(2) == 0);

  //---- Increase overlap along x-axis -----------------------------------------
    X_B.T() = EigenVector3<T>(2, 1, 0);

    succes = convex::growth_distance<T>(
        X_A, &A, X_B, &B, p_A, p_B, growth_scale, iterations, epsilon, max_its);
    BOOST_CHECK(succes);
    BOOST_CHECK_CLOSE(growth_scale, 0.5, 1.0);
    BOOST_CHECK_CLOSE(p_A(0), 1.5, 1.0);
    BOOST_CHECK_CLOSE(p_A(1), 1.0, 0.01);
    BOOST_CHECK(p_A(2) == 0);
    BOOST_CHECK_CLOSE(p_B(0), 1.5, 1.0);
    BOOST_CHECK_CLOSE(p_B(1), 1.0, 0.01);
    BOOST_CHECK(p_B(2) == 0);
}

BOOST_AUTO_TEST_CASE(growth_distance_capsule_and_sphere)
{
    convex::Capsule<T> A;
    geometry::Sphere<T> B;

    CoordSysEigen<T> X_A, X_B;
    EigenVector3<T> p_A, p_B;
    T growth_scale;
    size_t iterations;

  //--- Sphere touching capsule in one point in the middle of the --------------
  //--- "cylinder" face
    A.half_height() = 1.0;
    A.radius() = .50;
    B.radius() = .50;
    X_A.T() = EigenVector3<T>(.5, 0, 1);
    X_B.T() = EigenVector3<T>(1.5, 0, 1);

    T const epsilon = 0.01;

    bool succes = convex::growth_distance<T>(
        X_A, &A, X_B, &B, p_A, p_B, growth_scale, iterations, epsilon, 10u);

    BOOST_CHECK(succes);
    BOOST_CHECK_CLOSE(growth_scale, 1.0, 1.0);

    BOOST_CHECK_CLOSE(p_A(0), 1.0, 1.0);
    BOOST_CHECK_SMALL(p_A(1), 0.001);
    BOOST_CHECK_CLOSE(p_A(2), 1.0, 0.01);

    BOOST_CHECK_CLOSE(p_B(0), 1.0, 1.0);
    BOOST_CHECK_SMALL(p_B(1), 0.001);
    BOOST_CHECK_CLOSE(p_B(2), 1.0, 0.01);

  //--- sphere is overlaping top-shere cap of capsule --------------------------
    X_B.T() = EigenVector3<T>(.5, 0, 2);

    succes = convex::growth_distance<T>(X_A, &A, X_B, &B, p_A, p_B,
                                        growth_scale, iterations, epsilon, 10u);
    BOOST_CHECK(succes);
    BOOST_CHECK_CLOSE(growth_scale, 0.5, 1.0);

    BOOST_CHECK_CLOSE(p_A(0), 0.5, 0.01);
    BOOST_CHECK_CLOSE(p_A(1), 0.0, 0.01);
    BOOST_CHECK_CLOSE(p_A(2), 1.75, 1.0);

    BOOST_CHECK_CLOSE(p_B(0), 0.5, 0.01);
    BOOST_CHECK_CLOSE(p_B(1), 0.0, 0.01);
    BOOST_CHECK_CLOSE(p_B(2), 1.75, 1.0);

  //--- sphere is separated along y-axis at bottom-shere cap of capsule --------
    X_B.T() = EigenVector3<T>(.5, 2, 0);

    succes = convex::growth_distance<T>(X_A, &A, X_B, &B, p_A, p_B,
                                        growth_scale, iterations, 0.01, 10u);
    BOOST_CHECK(succes);
    BOOST_CHECK_CLOSE(growth_scale, 2.0, 1.0);
    BOOST_CHECK_CLOSE(p_A(0), 0.5, 0.01);
    BOOST_CHECK_CLOSE(p_A(1), 1.0, 1.0);
    BOOST_CHECK_SMALL(p_A(2), 0.01);

    BOOST_CHECK_CLOSE(p_B(0), 0.5, 0.01);
    BOOST_CHECK_CLOSE(p_B(1), 1.0, 1.0);
    BOOST_CHECK_SMALL(p_B(2), 0.01);
}

BOOST_AUTO_TEST_CASE(growth_distance_cylinders)
{
    convex::Cylinder<T> A;
    convex::Cylinder<T> B;

    CoordSysEigen<T> X_A, X_B;
    EigenVector3<T> p_A, p_B;
    T growth_scale;
    size_t iterations;

    A.half_height() = 2.0f;
    A.radius() = 1.0f;
    B.half_height() = 1.0f;
    B.radius() = 1.0f;

    X_A.T() = EigenVector3<T>(1, 0, 2);
    X_B.T() = EigenVector3<T>(2, 0, 2);

    //--- Cylinder's overlapping along-x-axis
    bool succes = convex::growth_distance<T>(
        X_A, &A, X_B, &B, p_A, p_B, growth_scale, iterations, 0.01f, 10u);

    BOOST_CHECK(succes);
    BOOST_CHECK_CLOSE(growth_scale, 0.5, 1.0f);

    BOOST_CHECK_CLOSE(p_A(0), 1.5, 1.0);
    BOOST_CHECK_SMALL(p_A(1), 0.001);
    BOOST_CHECK_CLOSE(p_A(2), 2.0, 0.01);

    BOOST_CHECK_CLOSE(p_B(0), 1.5, 1.0);
    BOOST_CHECK_SMALL(p_B(1), 0.01);
    BOOST_CHECK_CLOSE(p_B(2), 2.0, 0.01);

  //--- Cylinders on-top of each other ----------------------------------------
    X_B.T() = EigenVector3<T>(1, 0, 4);

    succes = convex::growth_distance<T>(X_A, &A, X_B, &B, p_A, p_B,
                                        growth_scale, iterations, 0.01f, 10u);
    BOOST_CHECK(succes);
    BOOST_CHECK_CLOSE(growth_scale, 0.6666666, 1.0f);

    BOOST_CHECK_CLOSE(p_A(0), 1.0, 1.0);
    BOOST_CHECK_SMALL(p_A(1), 0.001);
    BOOST_CHECK_CLOSE(p_A(2), 3.3333333, 1.0);

    BOOST_CHECK_CLOSE(p_B(0), 1.0, 1.0);
    BOOST_CHECK_SMALL(p_B(1), 0.01);
    BOOST_CHECK_CLOSE(p_B(2), 3.3333333, 1.0);

  //--- overlapping along x-y direction ----------------------------------------

    B.radius() = 2.0f;
    B.half_height() = 2.0f;
    X_B.T() = EigenVector3<T>(3, 0, 4);

    succes = convex::growth_distance<T>(X_A, &A, X_B, &B, p_A, p_B,
                                        growth_scale, iterations, 0.01f, 10u);

    BOOST_CHECK(succes);
    BOOST_CHECK_CLOSE(growth_scale, 0.6666666, 1.0f);

    BOOST_CHECK_CLOSE(p_A(0), 1.666666, 1.0);
    BOOST_CHECK_SMALL(p_A(1), 0.001);
    BOOST_CHECK_CLOSE(p_A(2), 3.0, 0.01);

    BOOST_CHECK_CLOSE(p_B(0), 1.666666, 1.0);
    BOOST_CHECK_SMALL(p_B(1), 0.01);
    BOOST_CHECK_CLOSE(p_B(2), 3.0, 0.01);
}

BOOST_AUTO_TEST_CASE(growth_distance_box_sphere)
{
    geometry::BoxEigen<T> A;
    geometry::Sphere<T> B;

    CoordSysEigen<T> X_A, X_B;
    EigenVector3<T> p_A = EigenVector3<T>(0, 0, 0);
    EigenVector3<T> p_B = EigenVector3<T>(0, 0, 0);
    T growth_scale;

    size_t iterations;
    size_t const max_its = 1000u;
    T const epsilon = 0.00001;

    A.half_extent() = EigenVector3<T>(.5, .5, .5);
    B.radius() = .5;

    X_A.T() = EigenVector3<T>(0, 0, 0);
    X_A.Q() = EigenQuaternion<T>::Identity();

    X_B.T() = EigenVector3<T>(0.0, 1.0, 0.0);
    X_B.Q() = EigenQuaternion<T>::Identity();

    bool succes = convex::growth_distance<T>(
        X_A, &A, X_B, &B, p_A, p_B, growth_scale, iterations, epsilon, max_its);

    BOOST_CHECK(succes);

    BOOST_CHECK_CLOSE(growth_scale, 1.0, 6.0f);

    BOOST_CHECK_SMALL(p_A(0), 0.05);
    BOOST_CHECK_CLOSE(p_A(1), 0.5, 6.0);
    BOOST_CHECK_SMALL(p_A(2), 0.05);

    BOOST_CHECK_SMALL(p_B(0), 0.05);
    BOOST_CHECK_CLOSE(p_B(1), 0.5, 6.0);
    BOOST_CHECK_SMALL(p_B(2), 0.05);

  //--- Sphere and box are separated by large distance ------------------------
    X_B.T() = EigenVector3<T>(1, 0, 4);

    succes = convex::growth_distance<T>(X_A, &A, X_B, &B, p_A, p_B,
                                        growth_scale, iterations, 0.01f, 10u);

    BOOST_CHECK(succes);
    BOOST_CHECK_CLOSE(growth_scale, 4.0, 5.0);

    BOOST_CHECK_CLOSE(p_A(0), 1.0,
                      25.0);   // 2015-12-11 Kenny: Really crapy results due to bad accuracy in inner GJK calls
    BOOST_CHECK_SMALL(p_A(1), 0.05);
    BOOST_CHECK_CLOSE(p_A(2), 2.0, 5.0);

    BOOST_CHECK_CLOSE(p_B(0), 1.0,
                      25.0);   // 2015-12-11 Kenny: Really crapy results due to bad accuracy in inner GJK calls
    BOOST_CHECK_SMALL(p_B(1), 0.05);
    BOOST_CHECK_CLOSE(p_B(2), 2.0, 5.0);
}

BOOST_AUTO_TEST_CASE(growth_distance_box_box_compile_test)
{
    geometry::BoxEigen<T> A;
    geometry::BoxEigen<T> B;

    CoordSysEigen<T> X_A, X_B;
    EigenVector3<T> p_A, p_B;
    T growth_scale;
    size_t iterations;

    A.half_extent() = EigenVector3<T>(1.0, 1.0f, 1.0f);
    B.half_extent() = EigenVector3<T>(1.0f, 1.0f, 1.0f);

    T delta = (0.01);
    T epsilon = std::min(delta * A.get_scale(), delta * B.get_scale());
    X_A.Q() = EigenQuaternion<T>(1, 0, 0, 0);
    X_B.Q() = EigenQuaternion<T>(1, 0, 0, 0);

    X_A.T() = EigenVector3<T>(0, -1.0, 0);
    X_B.T() = EigenVector3<T>(0, 1, 0);

    convex::growth_distance<T>(X_A, &A, X_B, &B, p_A, p_B, growth_scale,
                               iterations, epsilon, 100u);

    X_A.T() = EigenVector3<T>(0, -1.0, 0);
    X_B.T() = EigenVector3<T>(0, 2, 0);

    convex::growth_distance<T>(X_A, &A, X_B, &B, p_A, p_B, growth_scale,
                               iterations, epsilon, 100u);

    X_A.T() = EigenVector3<T>(0, -1.0, 0);
    X_B.T() = EigenVector3<T>(0, 0.9, 0);

    convex::growth_distance<T>(X_A, &A, X_B, &B, p_A, p_B, growth_scale,
                               iterations, epsilon, 100u);
}

BOOST_AUTO_TEST_SUITE_END();
