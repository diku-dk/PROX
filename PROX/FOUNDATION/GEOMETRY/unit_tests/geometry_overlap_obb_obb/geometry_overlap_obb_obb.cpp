#include <geometry.h>
#include <tiny.h>

#define BOOST_AUTO_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include <boost/test/unit_test_suite.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/test/test_tools.hpp>

using T = float;

BOOST_AUTO_TEST_SUITE(geometry);

BOOST_AUTO_TEST_CASE(overlap_obb_obb_test)
{
  // touching bottom-top faces
    {
        EigenVector3<T> const centerA = EigenVector3<T>(0.0, 0.0, 0.0);
        EigenVector3<T> const half_extA = EigenVector3<T>(1.0, 1.0, 1.0);
        EigenQuaternion<T> const qA = EigenQuaternion<T>::Identity();

        geometry::OBBEigen<T> obbA
            = geometry::make_obb<T>(centerA, qA, half_extA);

        EigenVector3<T> const centerB = EigenVector3<T>(0.0, -2.0, 0.0);
        EigenVector3<T> const half_extB = EigenVector3<T>(1.0, 1.0, 1.0);
        EigenQuaternion<T> const qB = Rotatey(std::numbers::pi_v<T> * 0.25f);

        geometry::OBBEigen<T> obbB
            = geometry::make_obb<T>(centerB, qB, half_extB);

        bool const test1 = geometry::overlap_obb_obb(obbA, obbB);
        bool const test2 = geometry::overlap_obb_obb(obbB, obbA);

        BOOST_CHECK(test1);
        BOOST_CHECK(test2);
    }
  // touching left-right faces
    {
        EigenVector3<T> const centerA = EigenVector3<T>(0.0, 0.0, 0.0);
        EigenVector3<T> const half_extA = EigenVector3<T>(1.0, 1.0, 1.0);
        EigenQuaternion<T> const qA = EigenQuaternion<T>::Identity();

        geometry::OBBEigen<T> obbA
            = geometry::make_obb<T>(centerA, qA, half_extA);

        EigenVector3<T> const centerB = EigenVector3<T>(2.0, 0.0, 0.0);
        EigenVector3<T> const half_extB = EigenVector3<T>(1.0, 1.0, 1.0);
        EigenQuaternion<T> const qB = Rotatex(std::numbers::pi_v<T> * 0.25f);

        geometry::OBBEigen<T> obbB
            = geometry::make_obb<T>(centerB, qB, half_extB);

        bool const test1 = geometry::overlap_obb_obb(obbA, obbB);
        bool const test2 = geometry::overlap_obb_obb(obbB, obbA);

        BOOST_CHECK(test1);
        BOOST_CHECK(test2);
    }
  // touching front-back faces
    {
        EigenVector3<T> const centerA = EigenVector3<T>(0.0, 0.0, 0.0);
        EigenVector3<T> const half_extA = EigenVector3<T>(1.0, 1.0, 1.0);
        EigenQuaternion<T> const qA = EigenQuaternion<T>::Identity();

        geometry::OBBEigen<T> obbA
            = geometry::make_obb<T>(centerA, qA, half_extA);

        EigenVector3<T> const centerB = EigenVector3<T>(0.0, 0.0, 2.0);
        EigenVector3<T> const half_extB = EigenVector3<T>(1.0, 1.0, 1.0);
        EigenQuaternion<T> const qB = Rotatez(std::numbers::pi_v<T> * 0.25f);

        geometry::OBBEigen<T> obbB
            = geometry::make_obb<T>(centerB, qB, half_extB);

        bool const test1 = geometry::overlap_obb_obb(obbA, obbB);
        bool const test2 = geometry::overlap_obb_obb(obbB, obbA);

        BOOST_CHECK(test1);
        BOOST_CHECK(test2);
    }
  // separating bottom-top faces
    {
        EigenVector3<T> const centerA = EigenVector3<T>(0.0, 0.0, 0.0);
        EigenVector3<T> const half_extA = EigenVector3<T>(1.0, 1.0, 1.0);
        EigenQuaternion<T> const qA = EigenQuaternion<T>::Identity();

        geometry::OBBEigen<T> obbA
            = geometry::make_obb<T>(centerA, qA, half_extA);

        EigenVector3<T> const centerB = EigenVector3<T>(0.0, -2.1, 0.0);
        EigenVector3<T> const half_extB = EigenVector3<T>(1.0, 1.0, 1.0);
        EigenQuaternion<T> const qB = Rotatey(std::numbers::pi_v<T> * 0.25f);

        geometry::OBBEigen<T> obbB
            = geometry::make_obb<T>(centerB, qB, half_extB);

        bool const test1 = geometry::overlap_obb_obb(obbA, obbB);
        bool const test2 = geometry::overlap_obb_obb(obbB, obbA);

        BOOST_CHECK(!test1);
        BOOST_CHECK(!test2);
    }
  // separating left-right faces
    {
        EigenVector3<T> const centerA = EigenVector3<T>(0.0, 0.0, 0.0);
        EigenVector3<T> const half_extA = EigenVector3<T>(1.0, 1.0, 1.0);
        EigenQuaternion<T> const qA = EigenQuaternion<T>::Identity();

        geometry::OBBEigen<T> obbA
            = geometry::make_obb<T>(centerA, qA, half_extA);

        EigenVector3<T> const centerB = EigenVector3<T>(2.1, 0.0, 0.0);
        EigenVector3<T> const half_extB = EigenVector3<T>(1.0, 1.0, 1.0);
        EigenQuaternion<T> const qB = Rotatex(std::numbers::pi_v<T> * 0.25f);

        geometry::OBBEigen<T> obbB
            = geometry::make_obb<T>(centerB, qB, half_extB);

        bool const test1 = geometry::overlap_obb_obb(obbA, obbB);
        bool const test2 = geometry::overlap_obb_obb(obbB, obbA);

        BOOST_CHECK(!test1);
        BOOST_CHECK(!test2);
    }
  // separating front-back faces
    {
        EigenVector3<T> const centerA = EigenVector3<T>(0.0, 0.0, 0.0);
        EigenVector3<T> const half_extA = EigenVector3<T>(1.0, 1.0, 1.0);
        EigenQuaternion<T> const qA = EigenQuaternion<T>::Identity();

        geometry::OBBEigen<T> obbA
            = geometry::make_obb<T>(centerA, qA, half_extA);

        EigenVector3<T> const centerB = EigenVector3<T>(0.0, 0.0, 2.1);
        EigenVector3<T> const half_extB = EigenVector3<T>(1.0, 1.0, 1.0);
        EigenQuaternion<T> const qB = Rotatez(std::numbers::pi_v<T> * 0.25f);

        geometry::OBBEigen<T> obbB
            = geometry::make_obb<T>(centerB, qB, half_extB);

        bool const test1 = geometry::overlap_obb_obb(obbA, obbB);
        bool const test2 = geometry::overlap_obb_obb(obbB, obbA);

        BOOST_CHECK(!test1);
        BOOST_CHECK(!test2);
    }
  // A inside of B
    {
        EigenVector3<T> const centerA = EigenVector3<T>(0.0, 0.0, 0.0);
        EigenVector3<T> const half_extA = EigenVector3<T>(0.5, 0.5, 0.5);
        EigenQuaternion<T> const qA = EigenQuaternion<T>::Identity();

        geometry::OBBEigen<T> obbA
            = geometry::make_obb<T>(centerA, qA, half_extA);

        EigenVector3<T> const centerB = EigenVector3<T>(0.0, 0.0, 0.0);
        EigenVector3<T> const half_extB = EigenVector3<T>(1.0, 1.0, 1.0);
        EigenQuaternion<T> const qB = Rotatey(std::numbers::pi_v<T> * 0.25f);

        geometry::OBBEigen<T> obbB
            = geometry::make_obb<T>(centerB, qB, half_extB);

        bool const test1 = geometry::overlap_obb_obb(obbA, obbB);
        bool const test2 = geometry::overlap_obb_obb(obbB, obbA);

        BOOST_CHECK(test1);
        BOOST_CHECK(test2);
    }
  // B inside of A
    {
        EigenVector3<T> const centerA = EigenVector3<T>(0.0, 0.0, 0.0);
        EigenVector3<T> const half_extA = EigenVector3<T>(1.0, 1.0, 1.0);
        EigenQuaternion<T> const qA = EigenQuaternion<T>::Identity();

        geometry::OBBEigen<T> obbA
            = geometry::make_obb<T>(centerA, qA, half_extA);

        EigenVector3<T> const centerB = EigenVector3<T>(0.0, 0.0, 0.0);
        EigenVector3<T> const half_extB = EigenVector3<T>(0.5, 0.5, 0.5);
        EigenQuaternion<T> const qB = Rotatey(std::numbers::pi_v<T> * 0.25f);

        geometry::OBBEigen<T> obbB
            = geometry::make_obb<T>(centerB, qB, half_extB);

        bool const test1 = geometry::overlap_obb_obb(obbA, obbB);
        bool const test2 = geometry::overlap_obb_obb(obbB, obbA);

        BOOST_CHECK(test1);
        BOOST_CHECK(test2);
    }
  // touching edge-edge-case
    {
        EigenVector3<T> const centerA = EigenVector3<T>(0.0, 0.0, 0.0);
        EigenVector3<T> const half_extA = EigenVector3<T>(1.0, 1.0, 1.0);
        EigenQuaternion<T> const qA = EigenQuaternion<T>::Identity();

        geometry::OBBEigen<T> obbA
            = geometry::make_obb<T>(centerA, qA, half_extA);

        EigenVector3<T> const centerB = EigenVector3<T>(2.0, 2.0, 0.0);
        EigenVector3<T> const half_extB = EigenVector3<T>(1.0, 1.0, 1.0);
        EigenQuaternion<T> const qB
            = Rotateu(std::numbers::pi_v<T> * 0.5f, centerB);

        geometry::OBBEigen<T> obbB
            = geometry::make_obb<T>(centerB, qB, half_extB);

        bool const test1 = geometry::overlap_obb_obb(obbA, obbB);
        bool const test2 = geometry::overlap_obb_obb(obbB, obbA);

        BOOST_CHECK(test1);
        BOOST_CHECK(test2);
    }
  // separating edge-edge-case
    {
        EigenVector3<T> const centerA = EigenVector3<T>(0.0, 0.0, 0.0);
        EigenVector3<T> const half_extA = EigenVector3<T>(1.0, 1.0, 1.0);
        EigenQuaternion<T> const qA = EigenQuaternion<T>::Identity();

        geometry::OBBEigen<T> obbA
            = geometry::make_obb<T>(centerA, qA, half_extA);

        EigenVector3<T> const centerB = EigenVector3<T>(2.01, 2.01, 0.0);
        EigenVector3<T> const half_extB = EigenVector3<T>(1.0, 1.0, 1.0);
        EigenQuaternion<T> const qB
            = Rotateu(std::numbers::pi_v<T> * 0.5f, centerB);

        geometry::OBBEigen<T> obbB
            = geometry::make_obb<T>(centerB, qB, half_extB);

        bool const test1 = geometry::overlap_obb_obb(obbA, obbB);
        bool const test2 = geometry::overlap_obb_obb(obbB, obbA);

        BOOST_CHECK(!test1);
        BOOST_CHECK(!test2);
    }
}

BOOST_AUTO_TEST_SUITE_END();
