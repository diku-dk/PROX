#include <geometry.h>
#include <tiny.h>

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

BOOST_AUTO_TEST_CASE(contacts_obb_sphere_test)
{
  // Touching rigth side
    {
        EigenVector3<T> const centerA = EigenVector3<T>(0.0, 0.0, 0.0);
        EigenVector3<T> const half_extA = EigenVector3<T>(1.0, 1.0, 1.0);
        EigenQuaternion<T> const qA = EigenQuaternion<T>::Identity();

        EigenVector3<T> const centerB = EigenVector3<T>(2.0, 0.0, 0.0);
        T const radiusB = 1.0;

        geometry::OBBEigen<T> const obb
            = geometry::make_obb<T>(centerA, qA, half_extA);
        geometry::Sphere<T> const sphere
            = geometry::make_sphere(centerB, radiusB);

        MyCallback callback;

        geometry::contacts_obb_sphere<T>(obb, sphere, 0.0, callback, false);

        BOOST_CHECK(callback.m_contacts.size() == 1u);

        BOOST_CHECK_CLOSE(callback.m_contacts[0].m_point[0], 1.0, 0.01);
        BOOST_CHECK_CLOSE(callback.m_contacts[0].m_point[1], 0.0, 0.01);
        BOOST_CHECK_CLOSE(callback.m_contacts[0].m_point[2], 0.0, 0.01);

        BOOST_CHECK_CLOSE(callback.m_contacts[0].m_normal[0], 1.0, 0.01);
        BOOST_CHECK_CLOSE(callback.m_contacts[0].m_normal[1], 0.0, 0.01);
        BOOST_CHECK_CLOSE(callback.m_contacts[0].m_normal[2], 0.0, 0.01);

        BOOST_CHECK_CLOSE(callback.m_contacts[0].m_distance, 0.0, 0.01);
    }
  // Separating rigth side
    {
        EigenVector3<T> const centerA = EigenVector3<T>(0.0, 0.0, 0.0);
        EigenVector3<T> const half_extA = EigenVector3<T>(1.0, 1.0, 1.0);
        EigenQuaternion<T> const qA = EigenQuaternion<T>::Identity();

        EigenVector3<T> const centerB = EigenVector3<T>(2.1, 0.0, 0.0);
        T const radiusB = 1.0;

        geometry::OBBEigen<T> const obb
            = geometry::make_obb<T>(centerA, qA, half_extA);
        geometry::Sphere<T> const sphere
            = geometry::make_sphere(centerB, radiusB);

        MyCallback callback;

        geometry::contacts_obb_sphere<T>(obb, sphere, 0.0, callback, false);

        BOOST_CHECK(callback.m_contacts.size() == 0u);
    }
  // Penetration rigth side
    {
        EigenVector3<T> const centerA = EigenVector3<T>(0.0, 0.0, 0.0);
        EigenVector3<T> const half_extA = EigenVector3<T>(1.0, 1.0, 1.0);
        EigenQuaternion<T> const qA = EigenQuaternion<T>::Identity();

        EigenVector3<T> const centerB = EigenVector3<T>(1.5, 0.0, 0.0);
        T const radiusB = 1.0;

        geometry::OBBEigen<T> const obb
            = geometry::make_obb<T>(centerA, qA, half_extA);
        geometry::Sphere<T> const sphere
            = geometry::make_sphere(centerB, radiusB);

        MyCallback callback;

        geometry::contacts_obb_sphere<T>(obb, sphere, 0.0, callback, false);

        BOOST_CHECK(callback.m_contacts.size() == 1u);

        BOOST_CHECK_CLOSE(callback.m_contacts[0].m_point[0], 1.0, 0.01);
        BOOST_CHECK_CLOSE(callback.m_contacts[0].m_point[1], 0.0, 0.01);
        BOOST_CHECK_CLOSE(callback.m_contacts[0].m_point[2], 0.0, 0.01);

        BOOST_CHECK_CLOSE(callback.m_contacts[0].m_normal[0], 1.0, 0.01);
        BOOST_CHECK_CLOSE(callback.m_contacts[0].m_normal[1], 0.0, 0.01);
        BOOST_CHECK_CLOSE(callback.m_contacts[0].m_normal[2], 0.0, 0.01);

        BOOST_CHECK_CLOSE(callback.m_contacts[0].m_distance, -0.5, 0.01);
    }
  // Penetration rigth side flipped-case
    {
        EigenVector3<T> const centerA = EigenVector3<T>(0.0, 0.0, 0.0);
        EigenVector3<T> const half_extA = EigenVector3<T>(1.0, 1.0, 1.0);
        EigenQuaternion<T> const qA = EigenQuaternion<T>::Identity();

        EigenVector3<T> const centerB = EigenVector3<T>(1.5, 0.0, 0.0);
        T const radiusB = 1.0;

        geometry::OBBEigen<T> const obb
            = geometry::make_obb<T>(centerA, qA, half_extA);
        geometry::Sphere<T> const sphere
            = geometry::make_sphere(centerB, radiusB);

        MyCallback callback;

        geometry::contacts_obb_sphere<T>(obb, sphere, 0.0, callback, true);

        BOOST_CHECK(callback.m_contacts.size() == 1u);

        BOOST_CHECK_CLOSE(callback.m_contacts[0].m_point[0], 1.0, 0.01);
        BOOST_CHECK_CLOSE(callback.m_contacts[0].m_point[1], 0.0, 0.01);
        BOOST_CHECK_CLOSE(callback.m_contacts[0].m_point[2], 0.0, 0.01);

        BOOST_CHECK_CLOSE(callback.m_contacts[0].m_normal[0], -1.0, 0.01);
        BOOST_CHECK_CLOSE(callback.m_contacts[0].m_normal[1], 0.0, 0.01);
        BOOST_CHECK_CLOSE(callback.m_contacts[0].m_normal[2], 0.0, 0.01);

        BOOST_CHECK_CLOSE(callback.m_contacts[0].m_distance, -0.5, 0.01);
    }

  // Touching left side
    {
        EigenVector3<T> const centerA = EigenVector3<T>(0.0, 0.0, 0.0);
        EigenVector3<T> const half_extA = EigenVector3<T>(1.0, 1.0, 1.0);
        EigenQuaternion<T> const qA = EigenQuaternion<T>::Identity();

        EigenVector3<T> const centerB = EigenVector3<T>(-2.0, 0.0, 0.0);
        T const radiusB = 1.0;

        geometry::OBBEigen<T> const obb
            = geometry::make_obb<T>(centerA, qA, half_extA);
        geometry::Sphere<T> const sphere
            = geometry::make_sphere(centerB, radiusB);

        MyCallback callback;

        geometry::contacts_obb_sphere<T>(obb, sphere, 0.0, callback, false);

        BOOST_CHECK(callback.m_contacts.size() == 1u);

        BOOST_CHECK_CLOSE(callback.m_contacts[0].m_point[0], -1.0, 0.01);
        BOOST_CHECK_CLOSE(callback.m_contacts[0].m_point[1], 0.0, 0.01);
        BOOST_CHECK_CLOSE(callback.m_contacts[0].m_point[2], 0.0, 0.01);

        BOOST_CHECK_CLOSE(callback.m_contacts[0].m_normal[0], -1.0, 0.01);
        BOOST_CHECK_CLOSE(callback.m_contacts[0].m_normal[1], 0.0, 0.01);
        BOOST_CHECK_CLOSE(callback.m_contacts[0].m_normal[2], 0.0, 0.01);

        BOOST_CHECK_CLOSE(callback.m_contacts[0].m_distance, 0.0, 0.01);
    }
  // Separating left side
    {
        EigenVector3<T> const centerA = EigenVector3<T>(0.0, 0.0, 0.0);
        EigenVector3<T> const half_extA = EigenVector3<T>(1.0, 1.0, 1.0);
        EigenQuaternion<T> const qA = EigenQuaternion<T>::Identity();

        EigenVector3<T> const centerB = EigenVector3<T>(-2.1, 0.0, 0.0);
        T const radiusB = 1.0;

        geometry::OBBEigen<T> const obb
            = geometry::make_obb<T>(centerA, qA, half_extA);
        geometry::Sphere<T> const sphere
            = geometry::make_sphere(centerB, radiusB);

        MyCallback callback;

        geometry::contacts_obb_sphere<T>(obb, sphere, 0.0, callback, false);

        BOOST_CHECK(callback.m_contacts.size() == 0u);
    }
  // Penetration left side
    {
        EigenVector3<T> const centerA = EigenVector3<T>(0.0, 0.0, 0.0);
        EigenVector3<T> const half_extA = EigenVector3<T>(1.0, 1.0, 1.0);
        EigenQuaternion<T> const qA = EigenQuaternion<T>::Identity();

        EigenVector3<T> const centerB = EigenVector3<T>(-1.5, 0.0, 0.0);
        T const radiusB = 1.0;

        geometry::OBBEigen<T> const obb
            = geometry::make_obb<T>(centerA, qA, half_extA);
        geometry::Sphere<T> const sphere
            = geometry::make_sphere(centerB, radiusB);

        MyCallback callback;

        geometry::contacts_obb_sphere<T>(obb, sphere, 0.0, callback, false);

        BOOST_CHECK(callback.m_contacts.size() == 1u);

        BOOST_CHECK_CLOSE(callback.m_contacts[0].m_point[0], -1.0, 0.01);
        BOOST_CHECK_CLOSE(callback.m_contacts[0].m_point[1], 0.0, 0.01);
        BOOST_CHECK_CLOSE(callback.m_contacts[0].m_point[2], 0.0, 0.01);

        BOOST_CHECK_CLOSE(callback.m_contacts[0].m_normal[0], -1.0, 0.01);
        BOOST_CHECK_CLOSE(callback.m_contacts[0].m_normal[1], 0.0, 0.01);
        BOOST_CHECK_CLOSE(callback.m_contacts[0].m_normal[2], 0.0, 0.01);

        BOOST_CHECK_CLOSE(callback.m_contacts[0].m_distance, -0.5, 0.01);
    }
  // Penetration left side flipped-case
    {
        EigenVector3<T> const centerA = EigenVector3<T>(0.0, 0.0, 0.0);
        EigenVector3<T> const half_extA = EigenVector3<T>(1.0, 1.0, 1.0);
        EigenQuaternion<T> const qA = EigenQuaternion<T>::Identity();

        EigenVector3<T> const centerB = EigenVector3<T>(-1.5, 0.0, 0.0);
        T const radiusB = 1.0;

        geometry::OBBEigen<T> const obb
            = geometry::make_obb<T>(centerA, qA, half_extA);
        geometry::Sphere<T> const sphere
            = geometry::make_sphere(centerB, radiusB);

        MyCallback callback;

        geometry::contacts_obb_sphere<T>(obb, sphere, 0.0, callback, true);

        BOOST_CHECK(callback.m_contacts.size() == 1u);

        BOOST_CHECK_CLOSE(callback.m_contacts[0].m_point[0], -1.0, 0.01);
        BOOST_CHECK_CLOSE(callback.m_contacts[0].m_point[1], 0.0, 0.01);
        BOOST_CHECK_CLOSE(callback.m_contacts[0].m_point[2], 0.0, 0.01);

        BOOST_CHECK_CLOSE(callback.m_contacts[0].m_normal[0], 1.0, 0.01);
        BOOST_CHECK_CLOSE(callback.m_contacts[0].m_normal[1], 0.0, 0.01);
        BOOST_CHECK_CLOSE(callback.m_contacts[0].m_normal[2], 0.0, 0.01);

        BOOST_CHECK_CLOSE(callback.m_contacts[0].m_distance, -0.5, 0.01);
    }
}

BOOST_AUTO_TEST_SUITE_END();
