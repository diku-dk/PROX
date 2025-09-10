#include <tiny.h>
#include <geometry.h> // needed for geometry::Sphere
#include <convex.h>

#define BOOST_AUTO_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include <boost/test/unit_test_suite.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/test/test_tools.hpp>

BOOST_AUTO_TEST_SUITE(motion_utils);

BOOST_AUTO_TEST_CASE(case_by_case_test)
{
/*    typedef tiny::MathTypes<double> M;
    typedef M::vector3_type V;
    typedef M::real_type T;
    typedef M::quaternion_type Q;
    typedef M::coordsys_type X;*/

    using T = double;
    CoordSysEigen<T> from;
    CoordSysEigen<T> to;

    from = CoordSysEigen<T>::identity();
    to.T() = EigenVector3<T>(1.0, 2.0, 3.0);
    to.Q() = Rotateu(3.0, EigenVector3<T>(1.0, 0.0, 0.0));
    EigenVector3<T> v;
    EigenVector3<T> w;
    convex::compute_velocities<T>(from, to, 1.0, v, w);

    BOOST_CHECK_CLOSE(v(0), to.T()(0), 0.01);
    BOOST_CHECK_CLOSE(v(1), to.T()(1), 0.01);
    BOOST_CHECK_CLOSE(v(2), to.T()(2), 0.01);
    BOOST_CHECK_CLOSE(w(0), 3.0, 0.01);
    BOOST_CHECK_CLOSE(w(1), 0.0, 0.01);
    BOOST_CHECK_CLOSE(w(2), 0.0, 0.01);

    T tau = 0.5;
    CoordSysEigen<T> cur = convex::integrate_motion<T>(from, tau, v, w);

    BOOST_CHECK_CLOSE(0.5, cur.T()(0), 0.01);
    BOOST_CHECK_CLOSE(1.0, cur.T()(1), 0.01);
    BOOST_CHECK_CLOSE(1.5, cur.T()(2), 0.01);

    T theta;
    EigenVector3<T> n;
    getAxisAngle(cur.Q(), n, theta);

    BOOST_CHECK_CLOSE(theta, 1.5, 0.01);
    BOOST_CHECK_CLOSE(n(0), 1.0, 0.01);
    BOOST_CHECK_CLOSE(n(1), 0.0, 0.01);
    BOOST_CHECK_CLOSE(n(2), 0.0, 0.01);

    geometry::Sphere<T> const A;
    geometry::Sphere<T> const B;
    CoordSysEigen<T> X_a;
    CoordSysEigen<T> X_b;

    X_a.T() = EigenVector3<T>(-2.0, 0.0, 0.0);
    X_a.Q() = EigenQuaternion<T>::Identity();

    X_b.T() = EigenVector3<T>(2.0, 0.0, 0.0);
    X_b.Q() = EigenQuaternion<T>::Identity();

    EigenVector3<T> p_a;
    EigenVector3<T> p_b;

    convex::compute_closest_points<T>(X_a, &A, X_b, &B, p_a, p_b);

    BOOST_CHECK_CLOSE(p_a(0), -1.0, 0.1);
    BOOST_CHECK_CLOSE(p_a(1), 0.0, 0.1);
    BOOST_CHECK_CLOSE(p_a(2), 0.0, 0.1);

    BOOST_CHECK_CLOSE(p_b(0), 1.0, 0.1);
    BOOST_CHECK_CLOSE(p_b(1), 0.0, 0.1);
    BOOST_CHECK_CLOSE(p_b(2), 0.0, 0.1);
}

BOOST_AUTO_TEST_SUITE_END();
