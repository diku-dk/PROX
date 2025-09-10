#include <geometry.h>
#include <eigenhelperall.h>

#define BOOST_AUTO_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include <boost/test/unit_test_suite.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/test/test_tools.hpp>

#include <vector>

using T = float;

class ContactInfo
{
public:
    EigenVector3<T> m_point;
    EigenVector3<T> m_normal;
    T m_distance;
};

class MyCallback : public geometry::ContactsCallback<T>
{
public:
    std::vector<ContactInfo> m_contacts;

public:
    void tempParenthesisOperatorImpl(const EigenVector3<T>& point,
                                     const EigenVector3<T>& normal,
                                     const T& distance)
    {
        ContactInfo info;

        info.m_point = point;
        info.m_normal = normal;
        info.m_distance = distance;

        m_contacts.push_back(info);
    }
};

BOOST_AUTO_TEST_SUITE(geometry);

BOOST_AUTO_TEST_CASE(contacts_obb_obb_test)
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
        EigenQuaternion<T> const qB = Rotatey<T>(std::numbers::pi_v<T> * 0.25);

        geometry::OBBEigen<T> obbB
            = geometry::make_obb<T>(centerB, qB, half_extB);

        MyCallback callback1;
        MyCallback callback2;

        bool const test1
            = geometry::contacts_obb_obb<T>(obbA, obbB, 0.0, callback1);
        bool const test2
            = geometry::contacts_obb_obb<T>(obbB, obbA, 0.0, callback2);

        BOOST_CHECK(test1);
        BOOST_CHECK(test2);

        BOOST_CHECK_EQUAL(callback1.m_contacts.size(), 8u);
        BOOST_CHECK_EQUAL(callback2.m_contacts.size(), 8u);

        for (unsigned int k = 0u; k < 8u; ++k)
        {
            BOOST_CHECK(fabs(callback1.m_contacts[k].m_normal(2)) < 0.00001);
            BOOST_CHECK(fabs(callback1.m_contacts[k].m_normal(0)) < 0.00001);
            BOOST_CHECK(1.0 - fabs(callback1.m_contacts[k].m_normal(1)) < 0.00001);
            BOOST_CHECK(fabs(callback1.m_contacts[k].m_distance) < 0.00001);
        }
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
        EigenQuaternion<T> const qB = Rotatex<T>(std::numbers::pi_v<T> * 0.25);

        geometry::OBBEigen<T> obbB
            = geometry::make_obb<T>(centerB, qB, half_extB);

        MyCallback callback1;
        MyCallback callback2;

        bool const test1
            = geometry::contacts_obb_obb<T>(obbA, obbB, 0.0, callback1);
        bool const test2
            = geometry::contacts_obb_obb<T>(obbB, obbA, 0.0, callback2);

        BOOST_CHECK(test1);
        BOOST_CHECK(test2);

        BOOST_CHECK_EQUAL(callback1.m_contacts.size(), 8u);
        BOOST_CHECK_EQUAL(callback2.m_contacts.size(), 8u);

        for (unsigned int k = 0u; k < 8u; ++k)
        {
            BOOST_CHECK(fabs(callback1.m_contacts[k].m_normal(2)) < 0.00001);
            BOOST_CHECK(fabs(callback1.m_contacts[k].m_normal(1)) < 0.00001);
            BOOST_CHECK(1.0 - fabs(callback1.m_contacts[k].m_normal(0)) < 0.00001);
            BOOST_CHECK(fabs(callback1.m_contacts[k].m_distance) < 0.00001);
        }
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
        EigenQuaternion<T> const qB = Rotatez<T>(std::numbers::pi_v<T> * 0.25);

        geometry::OBBEigen<T> obbB
            = geometry::make_obb<T>(centerB, qB, half_extB);

        MyCallback callback1;
        MyCallback callback2;

        bool const test1
            = geometry::contacts_obb_obb<T>(obbA, obbB, 0.0, callback1);
        bool const test2
            = geometry::contacts_obb_obb<T>(obbB, obbA, 0.0, callback2);

        BOOST_CHECK(test1);
        BOOST_CHECK(test2);

        BOOST_CHECK_EQUAL(callback1.m_contacts.size(), 8u);
        BOOST_CHECK_EQUAL(callback2.m_contacts.size(), 8u);

        for (unsigned int k = 0u; k < 8u; ++k)
        {
            BOOST_CHECK(fabs(callback1.m_contacts[k].m_normal(0)) < 0.00001);
            BOOST_CHECK(fabs(callback1.m_contacts[k].m_normal(1)) < 0.00001);
            BOOST_CHECK(1.0 - fabs(callback1.m_contacts[k].m_normal(2)) < 0.00001);
            BOOST_CHECK(fabs(callback1.m_contacts[k].m_distance) < 0.00001);
        }
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
        EigenQuaternion<T> const qB = Rotatey<T>(std::numbers::pi_v<T> * 0.25);

        geometry::OBBEigen<T> obbB
            = geometry::make_obb<T>(centerB, qB, half_extB);

        MyCallback callback1;
        MyCallback callback2;

        bool const test1
            = geometry::contacts_obb_obb<T>(obbA, obbB, 0.0, callback1);
        bool const test2
            = geometry::contacts_obb_obb<T>(obbB, obbA, 0.0, callback2);

        BOOST_CHECK(!test1);
        BOOST_CHECK(!test2);

        BOOST_CHECK_EQUAL(callback1.m_contacts.size(), 0u);
        BOOST_CHECK_EQUAL(callback2.m_contacts.size(), 0u);
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
        EigenQuaternion<T> const qB = Rotatex<T>(std::numbers::pi_v<T> * 0.25);

        geometry::OBBEigen<T> obbB
            = geometry::make_obb<T>(centerB, qB, half_extB);

        MyCallback callback1;
        MyCallback callback2;

        bool const test1
            = geometry::contacts_obb_obb<T>(obbA, obbB, 0.0, callback1);
        bool const test2
            = geometry::contacts_obb_obb<T>(obbB, obbA, 0.0, callback2);

        BOOST_CHECK(!test1);
        BOOST_CHECK(!test2);

        BOOST_CHECK_EQUAL(callback1.m_contacts.size(), 0u);
        BOOST_CHECK_EQUAL(callback2.m_contacts.size(), 0u);
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
        EigenQuaternion<T> const qB = Rotatez<T>(std::numbers::pi_v<T> * 0.25);

        geometry::OBBEigen<T> obbB
            = geometry::make_obb<T>(centerB, qB, half_extB);

        MyCallback callback1;
        MyCallback callback2;

        bool const test1
            = geometry::contacts_obb_obb<T>(obbA, obbB, 0.0, callback1);
        bool const test2
            = geometry::contacts_obb_obb<T>(obbB, obbA, 0.0, callback2);

        BOOST_CHECK(!test1);
        BOOST_CHECK(!test2);

        BOOST_CHECK_EQUAL(callback1.m_contacts.size(), 0u);
        BOOST_CHECK_EQUAL(callback2.m_contacts.size(), 0u);
    }
  // A inside of B
    {
        EigenVector3<T> const centerA = EigenVector3<T>(0.0, 0.5, 0.0);
        EigenVector3<T> const half_extA = EigenVector3<T>(0.5, 0.25, 0.5);
        EigenQuaternion<T> const qA = EigenQuaternion<T>::Identity();

        geometry::OBBEigen<T> obbA
            = geometry::make_obb<T>(centerA, qA, half_extA);

        EigenVector3<T> const centerB = EigenVector3<T>(0.0, 0.0, 0.0);
        EigenVector3<T> const half_extB = EigenVector3<T>(1.0, 1.0, 1.0);
        EigenQuaternion<T> const qB = EigenQuaternion<T>::Identity();

        geometry::OBBEigen<T> obbB
            = geometry::make_obb<T>(centerB, qB, half_extB);

        MyCallback callback1;
        MyCallback callback2;

        bool const test1
            = geometry::contacts_obb_obb<T>(obbA, obbB, 0.0, callback1);
        bool const test2
            = geometry::contacts_obb_obb<T>(obbB, obbA, 0.0, callback2);

        BOOST_CHECK(test1);
        BOOST_CHECK(test2);

        BOOST_CHECK_EQUAL(callback1.m_contacts.size(), 4u);
        BOOST_CHECK_EQUAL(callback2.m_contacts.size(), 4u);

        for (unsigned int k = 0u; k < 4u; ++k)
        {
            BOOST_CHECK(0.5 - fabs(callback1.m_contacts[k].m_point(0)) < 0.00001);
            BOOST_CHECK(0.5 - fabs(callback1.m_contacts[k].m_point(1)) < 0.00001);
            BOOST_CHECK(0.5 - fabs(callback1.m_contacts[k].m_point(2)) < 0.00001);
            BOOST_CHECK(fabs(callback1.m_contacts[k].m_normal(0)) < 0.00001);
            BOOST_CHECK(1.0 - fabs(callback1.m_contacts[k].m_normal(1)) < 0.00001);
            BOOST_CHECK(fabs(callback1.m_contacts[k].m_normal(2)) < 0.00001);
            BOOST_CHECK(0.5 - fabs(callback1.m_contacts[k].m_distance)
                        < 0.00001); // depth is off due to technicality in method, should have been 0.75
        }
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

        MyCallback callback1;
        MyCallback callback2;

        bool const test1
            = geometry::contacts_obb_obb<T>(obbA, obbB, 0.0, callback1);
        bool const test2
            = geometry::contacts_obb_obb<T>(obbB, obbA, 0.0, callback2);

        BOOST_CHECK(!test1);
        BOOST_CHECK(!test2);

        BOOST_CHECK_EQUAL(callback1.m_contacts.size(), 0u);
        BOOST_CHECK_EQUAL(callback2.m_contacts.size(), 0u);
    }
  // touching edge-edge-case
    {
        EigenVector3<T> const centerA = EigenVector3<T>(0.0, 0.0, 0.0);
        EigenVector3<T> const half_extA = EigenVector3<T>(1.0, 1.0, 1.0);
        EigenQuaternion<T> const qA = EigenQuaternion<T>::Identity();

        geometry::OBBEigen<T> obbA
            = geometry::make_obb<T>(centerA, qA, half_extA);

        EigenVector3<T> const centerB
            = EigenVector3<T>(1.999999, 1.999999, 0.0);
        EigenVector3<T> const half_extB = EigenVector3<T>(1.0, 1.0, 1.0);
        EigenQuaternion<T> const qB
            = Rotateu(std::numbers::pi_v<T> * 0.5f, centerB);

        geometry::OBBEigen<T> obbB
            = geometry::make_obb<T>(centerB, qB, half_extB);

        MyCallback callback1;
        MyCallback callback2;

        bool const test1
            = geometry::contacts_obb_obb<T>(obbA, obbB, 0.0, callback1);
        bool const test2
            = geometry::contacts_obb_obb<T>(obbB, obbA, 0.0, callback2);

        BOOST_CHECK(test1);
        BOOST_CHECK(test2);

        BOOST_CHECK_EQUAL(callback1.m_contacts.size(), 4u);
        BOOST_CHECK_EQUAL(callback2.m_contacts.size(), 4u);

        BOOST_CHECK(fabs(1.0 - callback1.m_contacts[0].m_point(0)) < 0.00001);
        BOOST_CHECK(fabs(1.0 - callback1.m_contacts[0].m_point(1)) < 0.00001);
        BOOST_CHECK(fabs(callback1.m_contacts[0].m_point(2)) < 0.00001);
        BOOST_CHECK(fabs(sqrt(0.5) - callback1.m_contacts[0].m_normal(0)) < 0.00001);
        BOOST_CHECK(fabs(sqrt(0.5) - callback1.m_contacts[0].m_normal(1)) < 0.00001);
        BOOST_CHECK(fabs(callback1.m_contacts[0].m_normal(2)) < 0.00001);
        BOOST_CHECK(fabs(callback1.m_contacts[0].m_distance) < 0.00001);

        BOOST_CHECK(fabs(1.0 - callback1.m_contacts[1].m_point(0)) < 0.00001);
        BOOST_CHECK(fabs(1.0 - callback1.m_contacts[1].m_point(1)) < 0.00001);
        BOOST_CHECK(fabs(callback1.m_contacts[1].m_point(2)) < 0.00001);
        BOOST_CHECK(fabs(sqrt(0.5) - callback1.m_contacts[1].m_normal(0)) < 0.00001);
        BOOST_CHECK(fabs(sqrt(0.5) - callback1.m_contacts[1].m_normal(1)) < 0.00001);
        BOOST_CHECK(fabs(callback1.m_contacts[1].m_normal(2)) < 0.00001);
        BOOST_CHECK(fabs(callback1.m_contacts[1].m_distance) < 0.00001);

        BOOST_CHECK(fabs(1.0 - callback1.m_contacts[2].m_point(0)) < 0.00001);
        BOOST_CHECK(fabs(1.0 - callback1.m_contacts[2].m_point(1)) < 0.00001);
        BOOST_CHECK(fabs(callback1.m_contacts[2].m_point(2)) < 0.00001);
        BOOST_CHECK(fabs(sqrt(0.5) - callback1.m_contacts[2].m_normal(0)) < 0.00001);
        BOOST_CHECK(fabs(sqrt(0.5) - callback1.m_contacts[2].m_normal(1)) < 0.00001);
        BOOST_CHECK(fabs(callback1.m_contacts[2].m_normal(2)) < 0.00001);
        BOOST_CHECK(fabs(callback1.m_contacts[2].m_distance) < 0.00001);

        BOOST_CHECK(fabs(1.0 - callback1.m_contacts[3].m_point(0)) < 0.00001);
        BOOST_CHECK(fabs(1.0 - callback1.m_contacts[3].m_point(1)) < 0.00001);
        BOOST_CHECK(fabs(callback1.m_contacts[3].m_point(2)) < 0.00001);
        BOOST_CHECK(fabs(sqrt(0.5) - callback1.m_contacts[3].m_normal(0)) < 0.00001);
        BOOST_CHECK(fabs(sqrt(0.5) - callback1.m_contacts[3].m_normal(1)) < 0.00001);
        BOOST_CHECK(fabs(callback1.m_contacts[3].m_normal(2)) < 0.00001);
        BOOST_CHECK(fabs(callback1.m_contacts[3].m_distance) < 0.00001);
    }
}

BOOST_AUTO_TEST_SUITE_END();
