#include <geometry.h>
#include <tiny.h>

#define BOOST_AUTO_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include <boost/test/unit_test_suite.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/test/test_tools.hpp>

using T = float;

BOOST_AUTO_TEST_SUITE(geometry);

BOOST_AUTO_TEST_CASE(obb_transform_test)
{
    EigenVector3<T> const center = EigenVector3<T>(1.0, 1.0, 1.0);
    EigenVector3<T> const half_ext = EigenVector3<T>(1.0, 2.0, 3.0);
    const EigenQuaternion<T> q = Rotatex(std::numbers::pi_v<T> * 0.5f);

    geometry::OBBEigen<T> const obb
        = geometry::make_obb<T>(center, q, half_ext);

    EigenVector3<T> const p0 = geometry::transform_to_obb(center, obb);

    BOOST_CHECK_CLOSE(p0(0), 0.0, 0.01);
    BOOST_CHECK_CLOSE(p0(1), 0.0, 0.01);
    BOOST_CHECK_CLOSE(p0(2), 0.0, 0.01);

    EigenVector3<T> const p1 = geometry::transform_from_obb(p0, obb);

    BOOST_CHECK_CLOSE(p1(0), 1.0, 0.01);
    BOOST_CHECK_CLOSE(p1(1), 1.0, 0.01);
    BOOST_CHECK_CLOSE(p1(2), 1.0, 0.01);

    EigenVector3<T> const p2
        = geometry::transform_from_obb(EigenVector3<T>(0.0, 2.0, 0.0), obb);

    BOOST_CHECK_CLOSE(p2(0), 1.0, 0.01);
    BOOST_CHECK_CLOSE(p2(1), 1.0, 0.01);
    BOOST_CHECK_CLOSE(p2(2), 3.0, 0.01);

    EigenVector3<T> const p3 = geometry::transform_to_obb(p2, obb);

    BOOST_CHECK_CLOSE(p3(0), 0.0, 0.01);
    BOOST_CHECK_CLOSE(p3(1), 2.0, 0.01);
    BOOST_CHECK_CLOSE(p3(2), 0.0, 0.01);
}

BOOST_AUTO_TEST_CASE(cylinder_transform_test)
{
    EigenVector3<T> const center = EigenVector3<T>(1.0, 1.0, 1.0);
    EigenVector3<T> const axis = EigenVector3<T>(1.0, 0.0, 0.0);
    T const radius = 1.0;
    T const height = 1.0;

    geometry::CylinderEigen<T> const cyl
        = geometry::make_cylinder(radius, height, axis, center);

    EigenVector3<T> const p0 = geometry::transform_to_cylinder(center, cyl);

    BOOST_CHECK_CLOSE(p0(0), 0.0, 0.01);
    BOOST_CHECK_CLOSE(p0(1), 0.0, 0.01);
    BOOST_CHECK_CLOSE(p0(2), 0.0, 0.01);

    EigenVector3<T> const p1 = geometry::transform_from_cylinder(p0, cyl);

    BOOST_CHECK_CLOSE(p1(0), 1.0, 0.01);
    BOOST_CHECK_CLOSE(p1(1), 1.0, 0.01);
    BOOST_CHECK_CLOSE(p1(2), 1.0, 0.01);

    EigenVector3<T> const p2 = geometry::transform_from_cylinder(
        EigenVector3<T>(0.0, 2.0, 0.0), cyl);

    BOOST_CHECK_CLOSE(p2(0), 1.0, 0.01);
    BOOST_CHECK_CLOSE(p2(1), 3.0, 0.01);
    BOOST_CHECK_CLOSE(p2(2), 1.0, 0.01);

    EigenVector3<T> const p3 = geometry::transform_to_cylinder(p2, cyl);

    BOOST_CHECK_CLOSE(p3(0), 0.0, 0.01);
    BOOST_CHECK_CLOSE(p3(1), 2.0, 0.01);
    BOOST_CHECK_CLOSE(p3(2), 0.0, 0.01);
}

BOOST_AUTO_TEST_SUITE_END();
