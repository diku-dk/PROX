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

BOOST_AUTO_TEST_CASE(contacts_sphere_tetrahedron_test)
{
  // Separation in vertex regions
    {
        EigenVector3<T> const center = EigenVector3<T>(2.1, 0.0, 0.0);
        T const radius = 1.0;

        EigenVector3<T> const p0 = EigenVector3<T>(0.0, 0.0, 0.0);
        EigenVector3<T> const p1 = EigenVector3<T>(1.0, 0.0, 0.0);
        EigenVector3<T> const p2 = EigenVector3<T>(0.0, 1.0, 0.0);
        EigenVector3<T> const p3 = EigenVector3<T>(0.0, 0.0, 1.0);

        geometry::Sphere<T> const sphere
            = geometry::make_sphere(center, radius);
        geometry::TetrahedronEigen<T> const tetrahedron
            = geometry::make_tetrahedron(p0, p1, p2, p3);

        MyCallback callback;

        bool const test = geometry::contacts_sphere_tetrahedron(sphere, tetrahedron, callback, false);

        BOOST_CHECK(!test);

        BOOST_CHECK_EQUAL(callback.m_contacts.size(), 0u);
    }
    {
        EigenVector3<T> const center = EigenVector3<T>(0.0, 2.1, 0.0);
        T const radius = 1.0;

        EigenVector3<T> const p0 = EigenVector3<T>(0.0, 0.0, 0.0);
        EigenVector3<T> const p1 = EigenVector3<T>(1.0, 0.0, 0.0);
        EigenVector3<T> const p2 = EigenVector3<T>(0.0, 1.0, 0.0);
        EigenVector3<T> const p3 = EigenVector3<T>(0.0, 0.0, 1.0);

        geometry::Sphere<T> const sphere
            = geometry::make_sphere(center, radius);
        geometry::TetrahedronEigen<T> const tetrahedron
            = geometry::make_tetrahedron(p0, p1, p2, p3);

        MyCallback callback;

        bool const test = geometry::contacts_sphere_tetrahedron(sphere, tetrahedron, callback, false);

        BOOST_CHECK(!test);

        BOOST_CHECK_EQUAL(callback.m_contacts.size(), 0u);
    }
    {
        EigenVector3<T> const center = EigenVector3<T>(0.0, 0.0, 2.1);
        T const radius = 1.0;

        EigenVector3<T> const p0 = EigenVector3<T>(0.0, 0.0, 0.0);
        EigenVector3<T> const p1 = EigenVector3<T>(1.0, 0.0, 0.0);
        EigenVector3<T> const p2 = EigenVector3<T>(0.0, 1.0, 0.0);
        EigenVector3<T> const p3 = EigenVector3<T>(0.0, 0.0, 1.0);

        geometry::Sphere<T> const sphere
            = geometry::make_sphere(center, radius);
        geometry::TetrahedronEigen<T> const tetrahedron
            = geometry::make_tetrahedron(p0, p1, p2, p3);

        MyCallback callback;

        bool const test = geometry::contacts_sphere_tetrahedron(sphere, tetrahedron, callback, false);

        BOOST_CHECK(!test);

        BOOST_CHECK_EQUAL(callback.m_contacts.size(), 0u);
    }
    {
        EigenVector3<T> const center = EigenVector3<T>(-1.0, -1.0, -1.0);
        T const radius = 1.0;

        EigenVector3<T> const p0 = EigenVector3<T>(0.0, 0.0, 0.0);
        EigenVector3<T> const p1 = EigenVector3<T>(1.0, 0.0, 0.0);
        EigenVector3<T> const p2 = EigenVector3<T>(0.0, 1.0, 0.0);
        EigenVector3<T> const p3 = EigenVector3<T>(0.0, 0.0, 1.0);

        geometry::Sphere<T> const sphere
            = geometry::make_sphere(center, radius);
        geometry::TetrahedronEigen<T> const tetrahedron
            = geometry::make_tetrahedron(p0, p1, p2, p3);

        MyCallback callback;

        bool const test = geometry::contacts_sphere_tetrahedron(sphere, tetrahedron, callback, false);

        BOOST_CHECK(!test);

        BOOST_CHECK_EQUAL(callback.m_contacts.size(), 0u);
    }
  // Separation in edge regions
    {
        EigenVector3<T> const center = EigenVector3<T>(0.5, -1.0, -0.5);
        T const radius = 1.0;

        EigenVector3<T> const p0 = EigenVector3<T>(0.0, 0.0, 0.0);
        EigenVector3<T> const p1 = EigenVector3<T>(1.0, 0.0, 0.0);
        EigenVector3<T> const p2 = EigenVector3<T>(0.0, 1.0, 0.0);
        EigenVector3<T> const p3 = EigenVector3<T>(0.0, 0.0, 1.0);

        geometry::Sphere<T> const sphere
            = geometry::make_sphere(center, radius);
        geometry::TetrahedronEigen<T> const tetrahedron
            = geometry::make_tetrahedron(p0, p1, p2, p3);

        MyCallback callback;

        bool const test = geometry::contacts_sphere_tetrahedron(sphere, tetrahedron, callback, false);

        BOOST_CHECK(!test);

        BOOST_CHECK_EQUAL(callback.m_contacts.size(), 0u);
    }
    {
        EigenVector3<T> const center = EigenVector3<T>(-1.0, 0.5, -0.5);
        T const radius = 1.0;

        EigenVector3<T> const p0 = EigenVector3<T>(0.0, 0.0, 0.0);
        EigenVector3<T> const p1 = EigenVector3<T>(1.0, 0.0, 0.0);
        EigenVector3<T> const p2 = EigenVector3<T>(0.0, 1.0, 0.0);
        EigenVector3<T> const p3 = EigenVector3<T>(0.0, 0.0, 1.0);

        geometry::Sphere<T> const sphere
            = geometry::make_sphere(center, radius);
        geometry::TetrahedronEigen<T> const tetrahedron
            = geometry::make_tetrahedron(p0, p1, p2, p3);

        MyCallback callback;

        bool const test = geometry::contacts_sphere_tetrahedron(sphere, tetrahedron, callback, false);

        BOOST_CHECK(!test);

        BOOST_CHECK_EQUAL(callback.m_contacts.size(), 0u);
    }
    {
        EigenVector3<T> const center = EigenVector3<T>(-1.0, -1.0, 0.5);
        T const radius = 1.0;

        EigenVector3<T> const p0 = EigenVector3<T>(0.0, 0.0, 0.0);
        EigenVector3<T> const p1 = EigenVector3<T>(1.0, 0.0, 0.0);
        EigenVector3<T> const p2 = EigenVector3<T>(0.0, 1.0, 0.0);
        EigenVector3<T> const p3 = EigenVector3<T>(0.0, 0.0, 1.0);

        geometry::Sphere<T> const sphere
            = geometry::make_sphere(center, radius);
        geometry::TetrahedronEigen<T> const tetrahedron
            = geometry::make_tetrahedron(p0, p1, p2, p3);

        MyCallback callback;

        bool const test = geometry::contacts_sphere_tetrahedron(sphere, tetrahedron, callback, false);

        BOOST_CHECK(!test);

        BOOST_CHECK_EQUAL(callback.m_contacts.size(), 0u);
    }
    {
        EigenVector3<T> const center = EigenVector3<T>(1.0, 1.0, -1.0);
        T const radius = 1.0;

        EigenVector3<T> const p0 = EigenVector3<T>(0.0, 0.0, 0.0);
        EigenVector3<T> const p1 = EigenVector3<T>(1.0, 0.0, 0.0);
        EigenVector3<T> const p2 = EigenVector3<T>(0.0, 1.0, 0.0);
        EigenVector3<T> const p3 = EigenVector3<T>(0.0, 0.0, 1.0);

        geometry::Sphere<T> const sphere
            = geometry::make_sphere(center, radius);
        geometry::TetrahedronEigen<T> const tetrahedron
            = geometry::make_tetrahedron(p0, p1, p2, p3);

        MyCallback callback;

        bool const test = geometry::contacts_sphere_tetrahedron(sphere, tetrahedron, callback, false);

        BOOST_CHECK(!test);

        BOOST_CHECK_EQUAL(callback.m_contacts.size(), 0u);
    }
    {
        EigenVector3<T> const center = EigenVector3<T>(1.0, -1.0, 1.0);
        T const radius = 1.0;

        EigenVector3<T> const p0 = EigenVector3<T>(0.0, 0.0, 0.0);
        EigenVector3<T> const p1 = EigenVector3<T>(1.0, 0.0, 0.0);
        EigenVector3<T> const p2 = EigenVector3<T>(0.0, 1.0, 0.0);
        EigenVector3<T> const p3 = EigenVector3<T>(0.0, 0.0, 1.0);

        geometry::Sphere<T> const sphere
            = geometry::make_sphere(center, radius);
        geometry::TetrahedronEigen<T> const tetrahedron
            = geometry::make_tetrahedron(p0, p1, p2, p3);

        MyCallback callback;

        bool const test = geometry::contacts_sphere_tetrahedron(sphere, tetrahedron, callback, false);

        BOOST_CHECK(!test);

        BOOST_CHECK_EQUAL(callback.m_contacts.size(), 0u);
    }
    {
        EigenVector3<T> const center = EigenVector3<T>(-1.0, 1.0, 1.0);
        T const radius = 1.0;

        EigenVector3<T> const p0 = EigenVector3<T>(0.0, 0.0, 0.0);
        EigenVector3<T> const p1 = EigenVector3<T>(1.0, 0.0, 0.0);
        EigenVector3<T> const p2 = EigenVector3<T>(0.0, 1.0, 0.0);
        EigenVector3<T> const p3 = EigenVector3<T>(0.0, 0.0, 1.0);

        geometry::Sphere<T> const sphere
            = geometry::make_sphere(center, radius);
        geometry::TetrahedronEigen<T> const tetrahedron
            = geometry::make_tetrahedron(p0, p1, p2, p3);

        MyCallback callback;

        bool const test = geometry::contacts_sphere_tetrahedron(sphere, tetrahedron, callback, false);

        BOOST_CHECK(!test);

        BOOST_CHECK_EQUAL(callback.m_contacts.size(), 0u);
    }
  // Separation in face regions
    {
        EigenVector3<T> const center = EigenVector3<T>(0.2, 0.2, -2.1);
        T const radius = 1.0;

        EigenVector3<T> const p0 = EigenVector3<T>(0.0, 0.0, 0.0);
        EigenVector3<T> const p1 = EigenVector3<T>(1.0, 0.0, 0.0);
        EigenVector3<T> const p2 = EigenVector3<T>(0.0, 1.0, 0.0);
        EigenVector3<T> const p3 = EigenVector3<T>(0.0, 0.0, 1.0);

        geometry::Sphere<T> const sphere
            = geometry::make_sphere(center, radius);
        geometry::TetrahedronEigen<T> const tetrahedron
            = geometry::make_tetrahedron(p0, p1, p2, p3);

        MyCallback callback;

        bool const test = geometry::contacts_sphere_tetrahedron(sphere, tetrahedron, callback, false);

        BOOST_CHECK(!test);

        BOOST_CHECK_EQUAL(callback.m_contacts.size(), 0u);
    }
    {
        EigenVector3<T> const center = EigenVector3<T>(0.2, -2.1, 0.2);
        T const radius = 1.0;

        EigenVector3<T> const p0 = EigenVector3<T>(0.0, 0.0, 0.0);
        EigenVector3<T> const p1 = EigenVector3<T>(1.0, 0.0, 0.0);
        EigenVector3<T> const p2 = EigenVector3<T>(0.0, 1.0, 0.0);
        EigenVector3<T> const p3 = EigenVector3<T>(0.0, 0.0, 1.0);

        geometry::Sphere<T> const sphere
            = geometry::make_sphere(center, radius);
        geometry::TetrahedronEigen<T> const tetrahedron
            = geometry::make_tetrahedron(p0, p1, p2, p3);

        MyCallback callback;

        bool const test = geometry::contacts_sphere_tetrahedron(sphere, tetrahedron, callback, false);

        BOOST_CHECK(!test);

        BOOST_CHECK_EQUAL(callback.m_contacts.size(), 0u);
    }
    {
        EigenVector3<T> const center = EigenVector3<T>(-2.1, 0.2, 0.2);
        T const radius = 1.0;

        EigenVector3<T> const p0 = EigenVector3<T>(0.0, 0.0, 0.0);
        EigenVector3<T> const p1 = EigenVector3<T>(1.0, 0.0, 0.0);
        EigenVector3<T> const p2 = EigenVector3<T>(0.0, 1.0, 0.0);
        EigenVector3<T> const p3 = EigenVector3<T>(0.0, 0.0, 1.0);

        geometry::Sphere<T> const sphere
            = geometry::make_sphere(center, radius);
        geometry::TetrahedronEigen<T> const tetrahedron
            = geometry::make_tetrahedron(p0, p1, p2, p3);

        MyCallback callback;

        bool const test = geometry::contacts_sphere_tetrahedron(sphere, tetrahedron, callback, false);

        BOOST_CHECK(!test);

        BOOST_CHECK_EQUAL(callback.m_contacts.size(), 0u);
    }
    {
        EigenVector3<T> const center = EigenVector3<T>(2.0, 2.0, 2.0);
        T const radius = 1.0;

        EigenVector3<T> const p0 = EigenVector3<T>(0.0, 0.0, 0.0);
        EigenVector3<T> const p1 = EigenVector3<T>(1.0, 0.0, 0.0);
        EigenVector3<T> const p2 = EigenVector3<T>(0.0, 1.0, 0.0);
        EigenVector3<T> const p3 = EigenVector3<T>(0.0, 0.0, 1.0);

        geometry::Sphere<T> const sphere
            = geometry::make_sphere(center, radius);
        geometry::TetrahedronEigen<T> const tetrahedron
            = geometry::make_tetrahedron(p0, p1, p2, p3);

        MyCallback callback;

        bool const test = geometry::contacts_sphere_tetrahedron(sphere, tetrahedron, callback, false);

        BOOST_CHECK(!test);

        BOOST_CHECK_EQUAL(callback.m_contacts.size(), 0u);
    }
  // Contacts in vertex regions
    {
        EigenVector3<T> const center = EigenVector3<T>(1.9, 0.0, 0.0);
        T const radius = 1.0;

        EigenVector3<T> const p0 = EigenVector3<T>(0.0, 0.0, 0.0);
        EigenVector3<T> const p1 = EigenVector3<T>(1.0, 0.0, 0.0);
        EigenVector3<T> const p2 = EigenVector3<T>(0.0, 1.0, 0.0);
        EigenVector3<T> const p3 = EigenVector3<T>(0.0, 0.0, 1.0);

        geometry::Sphere<T> const sphere
            = geometry::make_sphere(center, radius);
        geometry::TetrahedronEigen<T> const tetrahedron
            = geometry::make_tetrahedron(p0, p1, p2, p3);

        MyCallback callback;

        bool const test = geometry::contacts_sphere_tetrahedron(sphere, tetrahedron, callback, false);

        BOOST_CHECK(test);

        BOOST_CHECK_EQUAL(callback.m_contacts.size(), 1u);

        BOOST_CHECK_CLOSE(callback.m_contacts[0].m_point(0), 1.0, 0.01);
        BOOST_CHECK_CLOSE(callback.m_contacts[0].m_point(1), 0.0, 0.01);
        BOOST_CHECK_CLOSE(callback.m_contacts[0].m_point(2), 0.0, 0.01);

        BOOST_CHECK_CLOSE(callback.m_contacts[0].m_normal(0), -1.0, 0.01);
        BOOST_CHECK_CLOSE(callback.m_contacts[0].m_normal(1), 0.0, 0.01);
        BOOST_CHECK_CLOSE(callback.m_contacts[0].m_normal(2), 0.0, 0.01);

        BOOST_CHECK_CLOSE(callback.m_contacts[0].m_distance, -0.1, 0.01);
    }
    {
        EigenVector3<T> const center = EigenVector3<T>(1.9, 0.0, 0.0);
        T const radius = 1.0;

        EigenVector3<T> const p0 = EigenVector3<T>(0.0, 0.0, 0.0);
        EigenVector3<T> const p1 = EigenVector3<T>(1.0, 0.0, 0.0);
        EigenVector3<T> const p2 = EigenVector3<T>(0.0, 1.0, 0.0);
        EigenVector3<T> const p3 = EigenVector3<T>(0.0, 0.0, 1.0);

        geometry::Sphere<T> const sphere
            = geometry::make_sphere(center, radius);
        geometry::TetrahedronEigen<T> const tetrahedron
            = geometry::make_tetrahedron(p0, p1, p2, p3);

        MyCallback callback;

        bool const test = geometry::contacts_sphere_tetrahedron(sphere, tetrahedron, callback, true);

        BOOST_CHECK(test);

        BOOST_CHECK_EQUAL(callback.m_contacts.size(), 1u);

        BOOST_CHECK_CLOSE(callback.m_contacts[0].m_point(0), 1.0, 0.01);
        BOOST_CHECK_CLOSE(callback.m_contacts[0].m_point(1), 0.0, 0.01);
        BOOST_CHECK_CLOSE(callback.m_contacts[0].m_point(2), 0.0, 0.01);

        BOOST_CHECK_CLOSE(callback.m_contacts[0].m_normal(0), 1.0, 0.01);
        BOOST_CHECK_CLOSE(callback.m_contacts[0].m_normal(1), 0.0, 0.01);
        BOOST_CHECK_CLOSE(callback.m_contacts[0].m_normal(2), 0.0, 0.01);

        BOOST_CHECK_CLOSE(callback.m_contacts[0].m_distance, -0.1, 0.01);
    }

    {
        EigenVector3<T> const center = EigenVector3<T>(0.0, 1.9, 0.0);
        T const radius = 1.0;

        EigenVector3<T> const p0 = EigenVector3<T>(0.0, 0.0, 0.0);
        EigenVector3<T> const p1 = EigenVector3<T>(1.0, 0.0, 0.0);
        EigenVector3<T> const p2 = EigenVector3<T>(0.0, 1.0, 0.0);
        EigenVector3<T> const p3 = EigenVector3<T>(0.0, 0.0, 1.0);

        geometry::Sphere<T> const sphere
            = geometry::make_sphere(center, radius);
        geometry::TetrahedronEigen<T> const tetrahedron
            = geometry::make_tetrahedron(p0, p1, p2, p3);

        MyCallback callback;

        bool const test = geometry::contacts_sphere_tetrahedron(sphere, tetrahedron, callback, false);

        BOOST_CHECK(test);

        BOOST_CHECK_EQUAL(callback.m_contacts.size(), 1u);
    }
    {
        EigenVector3<T> const center = EigenVector3<T>(0.0, 0.0, 1.9);
        T const radius = 1.0;

        EigenVector3<T> const p0 = EigenVector3<T>(0.0, 0.0, 0.0);
        EigenVector3<T> const p1 = EigenVector3<T>(1.0, 0.0, 0.0);
        EigenVector3<T> const p2 = EigenVector3<T>(0.0, 1.0, 0.0);
        EigenVector3<T> const p3 = EigenVector3<T>(0.0, 0.0, 1.0);

        geometry::Sphere<T> const sphere
            = geometry::make_sphere(center, radius);
        geometry::TetrahedronEigen<T> const tetrahedron
            = geometry::make_tetrahedron(p0, p1, p2, p3);

        MyCallback callback;

        bool const test = geometry::contacts_sphere_tetrahedron(sphere, tetrahedron, callback, false);

        BOOST_CHECK(test);

        BOOST_CHECK_EQUAL(callback.m_contacts.size(), 1u);
    }
    {
        EigenVector3<T> const center = EigenVector3<T>(-0.5, -0.5, -0.5);
        T const radius = 1.0;

        EigenVector3<T> const p0 = EigenVector3<T>(0.0, 0.0, 0.0);
        EigenVector3<T> const p1 = EigenVector3<T>(1.0, 0.0, 0.0);
        EigenVector3<T> const p2 = EigenVector3<T>(0.0, 1.0, 0.0);
        EigenVector3<T> const p3 = EigenVector3<T>(0.0, 0.0, 1.0);

        geometry::Sphere<T> const sphere
            = geometry::make_sphere(center, radius);
        geometry::TetrahedronEigen<T> const tetrahedron
            = geometry::make_tetrahedron(p0, p1, p2, p3);

        MyCallback callback;

        bool const test = geometry::contacts_sphere_tetrahedron(sphere, tetrahedron, callback, false);

        BOOST_CHECK(test);

        BOOST_CHECK_EQUAL(callback.m_contacts.size(), 1u);
    }
  // Contact in edge regions
    {
        EigenVector3<T> const center = EigenVector3<T>(0.5, -0.5, -0.5);
        T const radius = 1.0;

        EigenVector3<T> const p0 = EigenVector3<T>(0.0, 0.0, 0.0);
        EigenVector3<T> const p1 = EigenVector3<T>(1.0, 0.0, 0.0);
        EigenVector3<T> const p2 = EigenVector3<T>(0.0, 1.0, 0.0);
        EigenVector3<T> const p3 = EigenVector3<T>(0.0, 0.0, 1.0);

        geometry::Sphere<T> const sphere
            = geometry::make_sphere(center, radius);
        geometry::TetrahedronEigen<T> const tetrahedron
            = geometry::make_tetrahedron(p0, p1, p2, p3);

        MyCallback callback;

        bool const test = geometry::contacts_sphere_tetrahedron(sphere, tetrahedron, callback, false);

        BOOST_CHECK(test);

        BOOST_CHECK_EQUAL(callback.m_contacts.size(), 1u);

        BOOST_CHECK_CLOSE(callback.m_contacts[0].m_point(0), 0.5, 0.01);
        BOOST_CHECK_CLOSE(callback.m_contacts[0].m_point(1), 0.0, 0.01);
        BOOST_CHECK_CLOSE(callback.m_contacts[0].m_point(2), 0.0, 0.01);

        BOOST_CHECK_CLOSE(callback.m_contacts[0].m_normal(0), 0.0, 0.01);
        BOOST_CHECK_CLOSE(callback.m_contacts[0].m_normal(1), 0.707106769, 0.01);
        BOOST_CHECK_CLOSE(callback.m_contacts[0].m_normal(2), 0.707106769, 0.01);

        BOOST_CHECK_CLOSE(callback.m_contacts[0].m_distance, -0.292893231, 0.01);
    }
    {
        EigenVector3<T> const center = EigenVector3<T>(0.5, -0.5, -0.5);
        T const radius = 1.0;

        EigenVector3<T> const p0 = EigenVector3<T>(0.0, 0.0, 0.0);
        EigenVector3<T> const p1 = EigenVector3<T>(1.0, 0.0, 0.0);
        EigenVector3<T> const p2 = EigenVector3<T>(0.0, 1.0, 0.0);
        EigenVector3<T> const p3 = EigenVector3<T>(0.0, 0.0, 1.0);

        geometry::Sphere<T> const sphere
            = geometry::make_sphere(center, radius);
        geometry::TetrahedronEigen<T> const tetrahedron
            = geometry::make_tetrahedron(p0, p1, p2, p3);

        MyCallback callback;

        bool const test = geometry::contacts_sphere_tetrahedron(sphere, tetrahedron, callback, true);

        BOOST_CHECK(test);

        BOOST_CHECK_EQUAL(callback.m_contacts.size(), 1u);

        BOOST_CHECK_CLOSE(callback.m_contacts[0].m_point(0), 0.5, 0.01);
        BOOST_CHECK_CLOSE(callback.m_contacts[0].m_point(1), 0.0, 0.01);
        BOOST_CHECK_CLOSE(callback.m_contacts[0].m_point(2), 0.0, 0.01);

        BOOST_CHECK_CLOSE(callback.m_contacts[0].m_normal(0), 0.0, 0.01);
        BOOST_CHECK_CLOSE(callback.m_contacts[0].m_normal(1), -0.707106769, 0.01);
        BOOST_CHECK_CLOSE(callback.m_contacts[0].m_normal(2), -0.707106769, 0.01);

        BOOST_CHECK_CLOSE(callback.m_contacts[0].m_distance, -0.292893231, 0.01);
    }
    {
        EigenVector3<T> const center = EigenVector3<T>(-0.5, 0.5, -0.5);
        T const radius = 1.0;

        EigenVector3<T> const p0 = EigenVector3<T>(0.0, 0.0, 0.0);
        EigenVector3<T> const p1 = EigenVector3<T>(1.0, 0.0, 0.0);
        EigenVector3<T> const p2 = EigenVector3<T>(0.0, 1.0, 0.0);
        EigenVector3<T> const p3 = EigenVector3<T>(0.0, 0.0, 1.0);

        geometry::Sphere<T> const sphere
            = geometry::make_sphere(center, radius);
        geometry::TetrahedronEigen<T> const tetrahedron
            = geometry::make_tetrahedron(p0, p1, p2, p3);

        MyCallback callback;

        bool const test = geometry::contacts_sphere_tetrahedron(sphere, tetrahedron, callback, false);

        BOOST_CHECK(test);

        BOOST_CHECK_EQUAL(callback.m_contacts.size(), 1u);
    }
    {
        EigenVector3<T> const center = EigenVector3<T>(-0.5, -0.5, 0.5);
        T const radius = 1.0;

        EigenVector3<T> const p0 = EigenVector3<T>(0.0, 0.0, 0.0);
        EigenVector3<T> const p1 = EigenVector3<T>(1.0, 0.0, 0.0);
        EigenVector3<T> const p2 = EigenVector3<T>(0.0, 1.0, 0.0);
        EigenVector3<T> const p3 = EigenVector3<T>(0.0, 0.0, 1.0);

        geometry::Sphere<T> const sphere
            = geometry::make_sphere(center, radius);
        geometry::TetrahedronEigen<T> const tetrahedron
            = geometry::make_tetrahedron(p0, p1, p2, p3);

        MyCallback callback;

        bool const test = geometry::contacts_sphere_tetrahedron(sphere, tetrahedron, callback, false);

        BOOST_CHECK(test);

        BOOST_CHECK_EQUAL(callback.m_contacts.size(), 1u);
    }
    {
        EigenVector3<T> const center = EigenVector3<T>(0.5, 0.5, -0.5);
        T const radius = 1.0;

        EigenVector3<T> const p0 = EigenVector3<T>(0.0, 0.0, 0.0);
        EigenVector3<T> const p1 = EigenVector3<T>(1.0, 0.0, 0.0);
        EigenVector3<T> const p2 = EigenVector3<T>(0.0, 1.0, 0.0);
        EigenVector3<T> const p3 = EigenVector3<T>(0.0, 0.0, 1.0);

        geometry::Sphere<T> const sphere
            = geometry::make_sphere(center, radius);
        geometry::TetrahedronEigen<T> const tetrahedron
            = geometry::make_tetrahedron(p0, p1, p2, p3);

        MyCallback callback;

        bool const test = geometry::contacts_sphere_tetrahedron(sphere, tetrahedron, callback, false);

        BOOST_CHECK(test);

        BOOST_CHECK_EQUAL(callback.m_contacts.size(), 1u);
    }
    {
        EigenVector3<T> const center = EigenVector3<T>(0.5, -0.5, 0.5);
        T const radius = 1.0;

        EigenVector3<T> const p0 = EigenVector3<T>(0.0, 0.0, 0.0);
        EigenVector3<T> const p1 = EigenVector3<T>(1.0, 0.0, 0.0);
        EigenVector3<T> const p2 = EigenVector3<T>(0.0, 1.0, 0.0);
        EigenVector3<T> const p3 = EigenVector3<T>(0.0, 0.0, 1.0);

        geometry::Sphere<T> const sphere
            = geometry::make_sphere(center, radius);
        geometry::TetrahedronEigen<T> const tetrahedron
            = geometry::make_tetrahedron(p0, p1, p2, p3);

        MyCallback callback;

        bool const test = geometry::contacts_sphere_tetrahedron(sphere, tetrahedron, callback, false);

        BOOST_CHECK(test);

        BOOST_CHECK_EQUAL(callback.m_contacts.size(), 1u);
    }
    {
        EigenVector3<T> const center = EigenVector3<T>(-0.5, 0.5, 0.5);
        T const radius = 1.0;

        EigenVector3<T> const p0 = EigenVector3<T>(0.0, 0.0, 0.0);
        EigenVector3<T> const p1 = EigenVector3<T>(1.0, 0.0, 0.0);
        EigenVector3<T> const p2 = EigenVector3<T>(0.0, 1.0, 0.0);
        EigenVector3<T> const p3 = EigenVector3<T>(0.0, 0.0, 1.0);

        geometry::Sphere<T> const sphere
            = geometry::make_sphere(center, radius);
        geometry::TetrahedronEigen<T> const tetrahedron
            = geometry::make_tetrahedron(p0, p1, p2, p3);

        MyCallback callback;

        bool const test = geometry::contacts_sphere_tetrahedron(sphere, tetrahedron, callback, false);

        BOOST_CHECK(test);

        BOOST_CHECK_EQUAL(callback.m_contacts.size(), 1u);
    }
  // Contact in face regions
    {
        EigenVector3<T> const center = EigenVector3<T>(0.2, 0.2, -0.9);
        T const radius = 1.0;

        EigenVector3<T> const p0 = EigenVector3<T>(0.0, 0.0, 0.0);
        EigenVector3<T> const p1 = EigenVector3<T>(1.0, 0.0, 0.0);
        EigenVector3<T> const p2 = EigenVector3<T>(0.0, 1.0, 0.0);
        EigenVector3<T> const p3 = EigenVector3<T>(0.0, 0.0, 1.0);

        geometry::Sphere<T> const sphere
            = geometry::make_sphere(center, radius);
        geometry::TetrahedronEigen<T> const tetrahedron
            = geometry::make_tetrahedron(p0, p1, p2, p3);

        MyCallback callback;

        bool const test = geometry::contacts_sphere_tetrahedron(sphere, tetrahedron, callback, false);

        BOOST_CHECK(test);

        BOOST_CHECK_EQUAL(callback.m_contacts.size(), 1u);

        BOOST_CHECK_CLOSE(callback.m_contacts[0].m_point(0), 0.2, 0.01);
        BOOST_CHECK_CLOSE(callback.m_contacts[0].m_point(1), 0.2, 0.01);
        BOOST_CHECK_CLOSE(callback.m_contacts[0].m_point(2), 0.0, 0.01);

        BOOST_CHECK_CLOSE(callback.m_contacts[0].m_normal(0), 0.0, 0.01);
        BOOST_CHECK_CLOSE(callback.m_contacts[0].m_normal(1), 0.0, 0.01);
        BOOST_CHECK_CLOSE(callback.m_contacts[0].m_normal(2), 1.0, 0.01);

        BOOST_CHECK_CLOSE(callback.m_contacts[0].m_distance, -0.1, 0.01);
    }
    {
        EigenVector3<T> const center = EigenVector3<T>(0.2, 0.2, -0.9);
        T const radius = 1.0;

        EigenVector3<T> const p0 = EigenVector3<T>(0.0, 0.0, 0.0);
        EigenVector3<T> const p1 = EigenVector3<T>(1.0, 0.0, 0.0);
        EigenVector3<T> const p2 = EigenVector3<T>(0.0, 1.0, 0.0);
        EigenVector3<T> const p3 = EigenVector3<T>(0.0, 0.0, 1.0);

        geometry::Sphere<T> const sphere
            = geometry::make_sphere(center, radius);
        geometry::TetrahedronEigen<T> const tetrahedron
            = geometry::make_tetrahedron(p0, p1, p2, p3);

        MyCallback callback;

        bool const test = geometry::contacts_sphere_tetrahedron(sphere, tetrahedron, callback, true);

        BOOST_CHECK(test);

        BOOST_CHECK_EQUAL(callback.m_contacts.size(), 1u);

        BOOST_CHECK_CLOSE(callback.m_contacts[0].m_point(0), 0.2, 0.01);
        BOOST_CHECK_CLOSE(callback.m_contacts[0].m_point(1), 0.2, 0.01);
        BOOST_CHECK_CLOSE(callback.m_contacts[0].m_point(2), 0.0, 0.01);

        BOOST_CHECK_CLOSE(callback.m_contacts[0].m_normal(0), 0.0, 0.01);
        BOOST_CHECK_CLOSE(callback.m_contacts[0].m_normal(1), 0.0, 0.01);
        BOOST_CHECK_CLOSE(callback.m_contacts[0].m_normal(2), -1.0, 0.01);

        BOOST_CHECK_CLOSE(callback.m_contacts[0].m_distance, -0.1, 0.01);
    }
    {
        EigenVector3<T> const center = EigenVector3<T>(0.2, -0.9, 0.2);
        T const radius = 1.0;

        EigenVector3<T> const p0 = EigenVector3<T>(0.0, 0.0, 0.0);
        EigenVector3<T> const p1 = EigenVector3<T>(1.0, 0.0, 0.0);
        EigenVector3<T> const p2 = EigenVector3<T>(0.0, 1.0, 0.0);
        EigenVector3<T> const p3 = EigenVector3<T>(0.0, 0.0, 1.0);

        geometry::Sphere<T> const sphere
            = geometry::make_sphere(center, radius);
        geometry::TetrahedronEigen<T> const tetrahedron
            = geometry::make_tetrahedron(p0, p1, p2, p3);

        MyCallback callback;

        bool const test = geometry::contacts_sphere_tetrahedron(sphere, tetrahedron, callback, false);

        BOOST_CHECK(test);

        BOOST_CHECK_EQUAL(callback.m_contacts.size(), 1u);
    }
    {
        EigenVector3<T> const center = EigenVector3<T>(-0.9, 0.2, 0.2);
        T const radius = 1.0;

        EigenVector3<T> const p0 = EigenVector3<T>(0.0, 0.0, 0.0);
        EigenVector3<T> const p1 = EigenVector3<T>(1.0, 0.0, 0.0);
        EigenVector3<T> const p2 = EigenVector3<T>(0.0, 1.0, 0.0);
        EigenVector3<T> const p3 = EigenVector3<T>(0.0, 0.0, 1.0);

        geometry::Sphere<T> const sphere
            = geometry::make_sphere(center, radius);
        geometry::TetrahedronEigen<T> const tetrahedron
            = geometry::make_tetrahedron(p0, p1, p2, p3);

        MyCallback callback;

        bool const test = geometry::contacts_sphere_tetrahedron(sphere, tetrahedron, callback, false);

        BOOST_CHECK(test);

        BOOST_CHECK_EQUAL(callback.m_contacts.size(), 1u);
    }
    {
        EigenVector3<T> const center = EigenVector3<T>(0.7, 0.7, 0.7);
        T const radius = 1.0;

        EigenVector3<T> const p0 = EigenVector3<T>(0.0, 0.0, 0.0);
        EigenVector3<T> const p1 = EigenVector3<T>(1.0, 0.0, 0.0);
        EigenVector3<T> const p2 = EigenVector3<T>(0.0, 1.0, 0.0);
        EigenVector3<T> const p3 = EigenVector3<T>(0.0, 0.0, 1.0);

        geometry::Sphere<T> const sphere
            = geometry::make_sphere(center, radius);
        geometry::TetrahedronEigen<T> const tetrahedron
            = geometry::make_tetrahedron(p0, p1, p2, p3);

        MyCallback callback;

        bool const test = geometry::contacts_sphere_tetrahedron(sphere, tetrahedron, callback, false);

        BOOST_CHECK(test);

        BOOST_CHECK_EQUAL(callback.m_contacts.size(), 1u);
    }

  // sphere inside
    {
        EigenVector3<T> const center = EigenVector3<T>(0.3, 0.2, 0.1);
        T const radius = 1.0;

        EigenVector3<T> const p0 = EigenVector3<T>(0.0, 0.0, 0.0);
        EigenVector3<T> const p1 = EigenVector3<T>(1.0, 0.0, 0.0);
        EigenVector3<T> const p2 = EigenVector3<T>(0.0, 1.0, 0.0);
        EigenVector3<T> const p3 = EigenVector3<T>(0.0, 0.0, 1.0);

        geometry::Sphere<T> const sphere
            = geometry::make_sphere(center, radius);
        geometry::TetrahedronEigen<T> const tetrahedron
            = geometry::make_tetrahedron(p0, p1, p2, p3);

        MyCallback callback;

        bool const test = geometry::contacts_sphere_tetrahedron(sphere, tetrahedron, callback, false);

        BOOST_CHECK(test);

        BOOST_CHECK_EQUAL(callback.m_contacts.size(), 1u);
    }
}

BOOST_AUTO_TEST_SUITE_END();
