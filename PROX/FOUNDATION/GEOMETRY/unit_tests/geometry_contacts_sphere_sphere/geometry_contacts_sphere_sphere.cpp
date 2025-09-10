#include "tiny_math_types.h"
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

BOOST_AUTO_TEST_CASE(contacts_sphere_sphere_test)
{

  // Penetration
    {
        EigenVector3<T> const centerA = EigenVector3<T>(0.0, 0.0, 0.0);
        T const radiusA = 2.0;
        EigenVector3<T> const centerB = EigenVector3<T>(5.0, 0.0, 0.0);
        T const radiusB = 4.0;

        geometry::Sphere<T> const A = geometry::make_sphere(centerA, radiusA);
        geometry::Sphere<T> const B = geometry::make_sphere(centerB, radiusB);

        MyCallback callback;

        geometry::contacts_sphere_sphere<T>(A, B, 0.0, callback);

        BOOST_CHECK(callback.m_contacts.size() == 1u);

        BOOST_CHECK_CLOSE(callback.m_contacts[0].m_point[0], 1.6666, 0.01);
        BOOST_CHECK_CLOSE(callback.m_contacts[0].m_point[1], 0.0, 0.01);
        BOOST_CHECK_CLOSE(callback.m_contacts[0].m_point[2], 0.0, 0.01);

        BOOST_CHECK_CLOSE(callback.m_contacts[0].m_normal[0], 1.0, 0.01);
        BOOST_CHECK_CLOSE(callback.m_contacts[0].m_normal[1], 0.0, 0.01);
        BOOST_CHECK_CLOSE(callback.m_contacts[0].m_normal[2], 0.0, 0.01);

        BOOST_CHECK_CLOSE(callback.m_contacts[0].m_distance, -1.0, 0.01);
    }

  // Touching
    {
        EigenVector3<T> const centerA = EigenVector3<T>(0.0, 0.0, 0.0);
        T const radiusA = 2.0;
        EigenVector3<T> const centerB = EigenVector3<T>(6.0, 0.0, 0.0);
        T const radiusB = 4.0;

        geometry::Sphere<T> const A = geometry::make_sphere(centerA, radiusA);
        geometry::Sphere<T> const B = geometry::make_sphere(centerB, radiusB);

        MyCallback callback;

        geometry::contacts_sphere_sphere<T>(A, B, 0.0, callback);

        BOOST_CHECK(callback.m_contacts.size() == 1u);

        BOOST_CHECK_CLOSE(callback.m_contacts[0].m_point[0], 2.0, 0.01);
        BOOST_CHECK_CLOSE(callback.m_contacts[0].m_point[1], 0.0, 0.01);
        BOOST_CHECK_CLOSE(callback.m_contacts[0].m_point[2], 0.0, 0.01);

        BOOST_CHECK_CLOSE(callback.m_contacts[0].m_normal[0], 1.0, 0.01);
        BOOST_CHECK_CLOSE(callback.m_contacts[0].m_normal[1], 0.0, 0.01);
        BOOST_CHECK_CLOSE(callback.m_contacts[0].m_normal[2], 0.0, 0.01);

        BOOST_CHECK_CLOSE(callback.m_contacts[0].m_distance, 0.0, 0.01);
    }

  // Separation
    {
        EigenVector3<T> const centerA = EigenVector3<T>(0.0, 0.0, 0.0);
        T const radiusA = 2.0;
        EigenVector3<T> const centerB = EigenVector3<T>(7.0, 0.0, 0.0);
        T const radiusB = 4.0;

        geometry::Sphere<T> const A = geometry::make_sphere(centerA, radiusA);
        geometry::Sphere<T> const B = geometry::make_sphere(centerB, radiusB);

        MyCallback callback;

        geometry::contacts_sphere_sphere<T>(A, B, 0.0, callback);

        BOOST_CHECK(callback.m_contacts.size() == 0u);
    }
}

BOOST_AUTO_TEST_SUITE_END();
