#include <geometry.h>
#include <tiny.h>

#define BOOST_AUTO_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include <boost/test/unit_test_suite.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/test/test_tools.hpp>

using T = float;

BOOST_AUTO_TEST_SUITE(geometry);

BOOST_AUTO_TEST_CASE(overlap_obb_capsule_test)
{
  // Separated on plus x-axis of OBB
    {
        EigenVector3<T> const center = EigenVector3<T>(0.0, 0.0, 0.0);
        EigenVector3<T> const half_ext = EigenVector3<T>(1.0, 1.0, 1.0);
        EigenQuaternion<T> const q = EigenQuaternion<T>::Identity();

        EigenVector3<T> const point0 = EigenVector3<T>(2.1, 0.0, 0.0);
        EigenVector3<T> const point1 = EigenVector3<T>(3.1, 2.1, 2.1);
        T const radius = 1.0;

        geometry::OBBEigen<T> const& obb
            = geometry::make_obb<T>(center, q, half_ext);
        geometry::CapsuleEigen<T> const& cap
            = geometry::make_capsule(radius, point0, point1);

        bool const test = geometry::overlap_obb_capsule(obb, cap);

        BOOST_CHECK(!test);
    }
  // Overlap from plus x-side of OBB
    {
        EigenVector3<T> const center = EigenVector3<T>(0.0, 0.0, 0.0);
        EigenVector3<T> const half_ext = EigenVector3<T>(1.0, 1.0, 1.0);
        EigenQuaternion<T> const q = EigenQuaternion<T>::Identity();

        EigenVector3<T> const point0 = EigenVector3<T>(1.9, 0.0, 0.0);
        EigenVector3<T> const point1 = EigenVector3<T>(3.1, 2.1, 2.1);
        T const radius = 1.0;

        geometry::OBBEigen<T> const& obb
            = geometry::make_obb<T>(center, q, half_ext);
        geometry::CapsuleEigen<T> const& cap
            = geometry::make_capsule(radius, point0, point1);

        bool const test = geometry::overlap_obb_capsule(obb, cap);

        BOOST_CHECK(test);
    }
  // Separated on negative x-axis of OBB
    {
        EigenVector3<T> const center = EigenVector3<T>(0.0, 0.0, 0.0);
        EigenVector3<T> const half_ext = EigenVector3<T>(1.0, 1.0, 1.0);
        EigenQuaternion<T> const q = EigenQuaternion<T>::Identity();

        EigenVector3<T> const point0 = EigenVector3<T>(-2.1, 0.0, 0.0);
        EigenVector3<T> const point1 = EigenVector3<T>(-3.1, 2.1, 2.1);
        T const radius = 1.0;

        geometry::OBBEigen<T> const& obb
            = geometry::make_obb<T>(center, q, half_ext);
        geometry::CapsuleEigen<T> const& cap
            = geometry::make_capsule(radius, point0, point1);

        bool const test = geometry::overlap_obb_capsule(obb, cap);

        BOOST_CHECK(!test);
    }
  // Overlap from negative x-side of OBB
    {
        EigenVector3<T> const center = EigenVector3<T>(0.0, 0.0, 0.0);
        EigenVector3<T> const half_ext = EigenVector3<T>(1.0, 1.0, 1.0);
        EigenQuaternion<T> const q = EigenQuaternion<T>::Identity();

        EigenVector3<T> const point0 = EigenVector3<T>(-1.9, 0.0, 0.0);
        EigenVector3<T> const point1 = EigenVector3<T>(-3.1, 2.1, 2.1);
        T const radius = 1.0;

        geometry::OBBEigen<T> const& obb
            = geometry::make_obb<T>(center, q, half_ext);
        geometry::CapsuleEigen<T> const& cap
            = geometry::make_capsule(radius, point0, point1);

        bool const test = geometry::overlap_obb_capsule(obb, cap);

        BOOST_CHECK(test);
    }
  // Separated on plus y-axis of OBB
    {
        EigenVector3<T> const center = EigenVector3<T>(0.0, 0.0, 0.0);
        EigenVector3<T> const half_ext = EigenVector3<T>(1.0, 1.0, 1.0);
        EigenQuaternion<T> const q = EigenQuaternion<T>::Identity();

        EigenVector3<T> const point0 = EigenVector3<T>(0.0, 2.1, 0.0);
        EigenVector3<T> const point1 = EigenVector3<T>(2.1, 3.1, 2.1);
        T const radius = 1.0;

        geometry::OBBEigen<T> const& obb
            = geometry::make_obb<T>(center, q, half_ext);
        geometry::CapsuleEigen<T> const& cap
            = geometry::make_capsule(radius, point0, point1);

        bool const test = geometry::overlap_obb_capsule(obb, cap);

        BOOST_CHECK(!test);
    }
  // Overlap from plus y-side of OBB
    {
        EigenVector3<T> const center = EigenVector3<T>(0.0, 0.0, 0.0);
        EigenVector3<T> const half_ext = EigenVector3<T>(1.0, 1.0, 1.0);
        EigenQuaternion<T> const q = EigenQuaternion<T>::Identity();

        EigenVector3<T> const point0 = EigenVector3<T>(0.0, 1.9, 0.0);
        EigenVector3<T> const point1 = EigenVector3<T>(2.1, 3.1, 2.1);
        T const radius = 1.0;

        geometry::OBBEigen<T> const& obb
            = geometry::make_obb<T>(center, q, half_ext);
        geometry::CapsuleEigen<T> const& cap
            = geometry::make_capsule(radius, point0, point1);

        bool const test = geometry::overlap_obb_capsule(obb, cap);

        BOOST_CHECK(test);
    }
  // Separated on negative y-axis of OBB
    {
        EigenVector3<T> const center = EigenVector3<T>(0.0, 0.0, 0.0);
        EigenVector3<T> const half_ext = EigenVector3<T>(1.0, 1.0, 1.0);
        EigenQuaternion<T> const q = EigenQuaternion<T>::Identity();

        EigenVector3<T> const point0 = EigenVector3<T>(0.0, -2.1, 0.0);
        EigenVector3<T> const point1 = EigenVector3<T>(2.1, -3.1, 2.1);
        T const radius = 1.0;

        geometry::OBBEigen<T> const& obb
            = geometry::make_obb<T>(center, q, half_ext);
        geometry::CapsuleEigen<T> const& cap
            = geometry::make_capsule(radius, point0, point1);

        bool const test = geometry::overlap_obb_capsule(obb, cap);

        BOOST_CHECK(!test);
    }
  // Overlap from negative y-side of OBB
    {
        EigenVector3<T> const center = EigenVector3<T>(0.0, 0.0, 0.0);
        EigenVector3<T> const half_ext = EigenVector3<T>(1.0, 1.0, 1.0);
        EigenQuaternion<T> const q = EigenQuaternion<T>::Identity();

        EigenVector3<T> const point0 = EigenVector3<T>(0.0, -1.9, 0.0);
        EigenVector3<T> const point1 = EigenVector3<T>(2.1, -3.1, 2.1);
        T const radius = 1.0;

        geometry::OBBEigen<T> const& obb
            = geometry::make_obb<T>(center, q, half_ext);
        geometry::CapsuleEigen<T> const& cap
            = geometry::make_capsule(radius, point0, point1);

        bool const test = geometry::overlap_obb_capsule(obb, cap);

        BOOST_CHECK(test);
    }
  // Separated on plus z-axis of OBB
    {
        EigenVector3<T> const center = EigenVector3<T>(0.0, 0.0, 0.0);
        EigenVector3<T> const half_ext = EigenVector3<T>(1.0, 1.0, 1.0);
        EigenQuaternion<T> const q = EigenQuaternion<T>::Identity();

        EigenVector3<T> const point0 = EigenVector3<T>(0.0, 0.0, 2.1);
        EigenVector3<T> const point1 = EigenVector3<T>(2.1, 2.1, 3.1);
        T const radius = 1.0;

        geometry::OBBEigen<T> const& obb
            = geometry::make_obb<T>(center, q, half_ext);
        geometry::CapsuleEigen<T> const& cap
            = geometry::make_capsule(radius, point0, point1);

        bool const test = geometry::overlap_obb_capsule(obb, cap);

        BOOST_CHECK(!test);
    }
  // Overlap from plus z-side of OBB
    {
        EigenVector3<T> const center = EigenVector3<T>(0.0, 0.0, 0.0);
        EigenVector3<T> const half_ext = EigenVector3<T>(1.0, 1.0, 1.0);
        EigenQuaternion<T> const q = EigenQuaternion<T>::Identity();

        EigenVector3<T> const point0 = EigenVector3<T>(0.0, 0.0, 1.9);
        EigenVector3<T> const point1 = EigenVector3<T>(2.1, 2.1, 3.1);
        T const radius = 1.0;

        geometry::OBBEigen<T> const& obb
            = geometry::make_obb<T>(center, q, half_ext);
        geometry::CapsuleEigen<T> const& cap
            = geometry::make_capsule(radius, point0, point1);

        bool const test = geometry::overlap_obb_capsule(obb, cap);

        BOOST_CHECK(test);
    }
  // Separated on negative z-axis of OBB
    {
        EigenVector3<T> const center = EigenVector3<T>(0.0, 0.0, 0.0);
        EigenVector3<T> const half_ext = EigenVector3<T>(1.0, 1.0, 1.0);
        EigenQuaternion<T> const q = EigenQuaternion<T>::Identity();

        EigenVector3<T> const point0 = EigenVector3<T>(0.0, 0.0, -2.1);
        EigenVector3<T> const point1 = EigenVector3<T>(2.1, 2.1, -3.1);
        T const radius = 1.0;

        geometry::OBBEigen<T> const& obb
            = geometry::make_obb<T>(center, q, half_ext);
        geometry::CapsuleEigen<T> const& cap
            = geometry::make_capsule(radius, point0, point1);

        bool const test = geometry::overlap_obb_capsule(obb, cap);

        BOOST_CHECK(!test);
    }
  // Overlap from negative z-side of OBB
    {
        EigenVector3<T> const center = EigenVector3<T>(0.0, 0.0, 0.0);
        EigenVector3<T> const half_ext = EigenVector3<T>(1.0, 1.0, 1.0);
        EigenQuaternion<T> const q = EigenQuaternion<T>::Identity();

        EigenVector3<T> const point0 = EigenVector3<T>(0.0, 0.0, -1.9);
        EigenVector3<T> const point1 = EigenVector3<T>(2.1, 2.1, -3.1);
        T const radius = 1.0;

        geometry::OBBEigen<T> const& obb
            = geometry::make_obb<T>(center, q, half_ext);
        geometry::CapsuleEigen<T> const& cap
            = geometry::make_capsule(radius, point0, point1);

        bool const test = geometry::overlap_obb_capsule(obb, cap);

        BOOST_CHECK(test);
    }
  // Separated by z-OBB X capsule axis
    {
        EigenVector3<T> const center = EigenVector3<T>(0.0, 0.0, 0.0);
        EigenVector3<T> const half_ext = EigenVector3<T>(1.0, 1.0, 1.0);
        EigenQuaternion<T> const q = EigenQuaternion<T>::Identity();

        EigenVector3<T> const point0 = EigenVector3<T>(0.0, -3.0, 0.0);
        EigenVector3<T> const point1 = EigenVector3<T>(3.0, 0.0, 0.0);
        T const radius = 0.5;

        geometry::OBBEigen<T> const& obb
            = geometry::make_obb<T>(center, q, half_ext);
        geometry::CapsuleEigen<T> const& cap
            = geometry::make_capsule(radius, point0, point1);

        bool const test = geometry::overlap_obb_capsule(obb, cap);

        BOOST_CHECK(!test);
    }
  // Overlap by z-OBB X capsule axis
    {
        EigenVector3<T> const center = EigenVector3<T>(0.0, 0.0, 0.0);
        EigenVector3<T> const half_ext = EigenVector3<T>(1.0, 1.0, 1.0);
        EigenQuaternion<T> const q = EigenQuaternion<T>::Identity();

        EigenVector3<T> const point0 = EigenVector3<T>(0.0, -3.0, 0.0);
        EigenVector3<T> const point1 = EigenVector3<T>(3.0, 0.0, 0.0);
        T const radius = 1.0;

        geometry::OBBEigen<T> const& obb
            = geometry::make_obb<T>(center, q, half_ext);
        geometry::CapsuleEigen<T> const& cap
            = geometry::make_capsule(radius, point0, point1);

        bool const test = geometry::overlap_obb_capsule(obb, cap);

        BOOST_CHECK(test);
    }
}

BOOST_AUTO_TEST_SUITE_END();
