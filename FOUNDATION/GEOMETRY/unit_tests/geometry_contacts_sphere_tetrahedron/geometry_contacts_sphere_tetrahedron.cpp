#include <geometry.h>
#include <tiny.h>

#define BOOST_AUTO_TEST_MAIN
#include <boost/test/auto_unit_test.hpp>
#include <boost/test/unit_test_suite.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/test_tools.hpp>


#include <vector>

typedef tiny::MathTypes<float> MT;
typedef MT::vector3_type       V;
typedef MT::quaternion_type    Q;
typedef MT::real_type          T;

class ContactInfo
{
public:

  V m_point;
  V m_normal;
  T m_distance;

};


class MyCallback
  : public geometry::ContactsCallback<V>
{
public:

  std::vector<ContactInfo> m_contacts;

public:

  void operator()(
                  V const & point
                  , V const & normal
                  , typename V::real_type const & distance
                  )
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
    V const center = V::make(2.1, 0.0, 0.0);
    T const radius = 1.0;

    V const p0 = V::make(0.0, 0.0, 0.0);
    V const p1 = V::make(1.0, 0.0, 0.0);
    V const p2 = V::make(0.0, 1.0, 0.0);
    V const p3 = V::make(0.0, 0.0, 1.0);

    geometry::Sphere<V>      const sphere      = geometry::make_sphere(center, radius);
    geometry::Tetrahedron<V> const tetrahedron = geometry::make_tetrahedron(p0, p1, p2 ,p3);

    MyCallback callback;

    bool const test = geometry::contacts_sphere_tetrahedron(sphere, tetrahedron, callback, false);

    BOOST_CHECK(!test);

    BOOST_CHECK_EQUAL(callback.m_contacts.size(), 0u);
  }
  {
    V const center = V::make(0.0, 2.1, 0.0);
    T const radius = 1.0;

    V const p0 = V::make(0.0, 0.0, 0.0);
    V const p1 = V::make(1.0, 0.0, 0.0);
    V const p2 = V::make(0.0, 1.0, 0.0);
    V const p3 = V::make(0.0, 0.0, 1.0);

    geometry::Sphere<V>      const sphere      = geometry::make_sphere(center, radius);
    geometry::Tetrahedron<V> const tetrahedron = geometry::make_tetrahedron(p0, p1, p2 ,p3);

    MyCallback callback;

    bool const test = geometry::contacts_sphere_tetrahedron(sphere, tetrahedron, callback, false);

    BOOST_CHECK(!test);

    BOOST_CHECK_EQUAL(callback.m_contacts.size(), 0u);
  }
  {
    V const center = V::make(0.0, 0.0, 2.1);
    T const radius = 1.0;

    V const p0 = V::make(0.0, 0.0, 0.0);
    V const p1 = V::make(1.0, 0.0, 0.0);
    V const p2 = V::make(0.0, 1.0, 0.0);
    V const p3 = V::make(0.0, 0.0, 1.0);

    geometry::Sphere<V>      const sphere      = geometry::make_sphere(center, radius);
    geometry::Tetrahedron<V> const tetrahedron = geometry::make_tetrahedron(p0, p1, p2 ,p3);

    MyCallback callback;

    bool const test = geometry::contacts_sphere_tetrahedron(sphere, tetrahedron, callback, false);

    BOOST_CHECK(!test);

    BOOST_CHECK_EQUAL(callback.m_contacts.size(), 0u);
  }
  {
    V const center = V::make(-1.0, -1.0, -1.0);
    T const radius = 1.0;

    V const p0 = V::make(0.0, 0.0, 0.0);
    V const p1 = V::make(1.0, 0.0, 0.0);
    V const p2 = V::make(0.0, 1.0, 0.0);
    V const p3 = V::make(0.0, 0.0, 1.0);

    geometry::Sphere<V>      const sphere      = geometry::make_sphere(center, radius);
    geometry::Tetrahedron<V> const tetrahedron = geometry::make_tetrahedron(p0, p1, p2 ,p3);

    MyCallback callback;

    bool const test = geometry::contacts_sphere_tetrahedron(sphere, tetrahedron, callback, false);

    BOOST_CHECK(!test);

    BOOST_CHECK_EQUAL(callback.m_contacts.size(), 0u);
  }
  // Separation in edge regions
  {
    V const center = V::make(0.5, -1.0, -0.5);
    T const radius = 1.0;

    V const p0 = V::make(0.0, 0.0, 0.0);
    V const p1 = V::make(1.0, 0.0, 0.0);
    V const p2 = V::make(0.0, 1.0, 0.0);
    V const p3 = V::make(0.0, 0.0, 1.0);

    geometry::Sphere<V>      const sphere      = geometry::make_sphere(center, radius);
    geometry::Tetrahedron<V> const tetrahedron = geometry::make_tetrahedron(p0, p1, p2 ,p3);

    MyCallback callback;

    bool const test = geometry::contacts_sphere_tetrahedron(sphere, tetrahedron, callback, false);

    BOOST_CHECK(!test);

    BOOST_CHECK_EQUAL(callback.m_contacts.size(), 0u);
  }
  {
    V const center = V::make(-1.0, 0.5, -0.5);
    T const radius = 1.0;

    V const p0 = V::make(0.0, 0.0, 0.0);
    V const p1 = V::make(1.0, 0.0, 0.0);
    V const p2 = V::make(0.0, 1.0, 0.0);
    V const p3 = V::make(0.0, 0.0, 1.0);

    geometry::Sphere<V>      const sphere      = geometry::make_sphere(center, radius);
    geometry::Tetrahedron<V> const tetrahedron = geometry::make_tetrahedron(p0, p1, p2 ,p3);

    MyCallback callback;

    bool const test = geometry::contacts_sphere_tetrahedron(sphere, tetrahedron, callback, false);

    BOOST_CHECK(!test);

    BOOST_CHECK_EQUAL(callback.m_contacts.size(), 0u);
  }
  {
    V const center = V::make(-1.0, -1.0, 0.5);
    T const radius = 1.0;

    V const p0 = V::make(0.0, 0.0, 0.0);
    V const p1 = V::make(1.0, 0.0, 0.0);
    V const p2 = V::make(0.0, 1.0, 0.0);
    V const p3 = V::make(0.0, 0.0, 1.0);

    geometry::Sphere<V>      const sphere      = geometry::make_sphere(center, radius);
    geometry::Tetrahedron<V> const tetrahedron = geometry::make_tetrahedron(p0, p1, p2 ,p3);

    MyCallback callback;

    bool const test = geometry::contacts_sphere_tetrahedron(sphere, tetrahedron, callback, false);

    BOOST_CHECK(!test);

    BOOST_CHECK_EQUAL(callback.m_contacts.size(), 0u);
  }
  {
    V const center = V::make( 1.0, 1.0, -1.0);
    T const radius = 1.0;

    V const p0 = V::make(0.0, 0.0, 0.0);
    V const p1 = V::make(1.0, 0.0, 0.0);
    V const p2 = V::make(0.0, 1.0, 0.0);
    V const p3 = V::make(0.0, 0.0, 1.0);

    geometry::Sphere<V>      const sphere      = geometry::make_sphere(center, radius);
    geometry::Tetrahedron<V> const tetrahedron = geometry::make_tetrahedron(p0, p1, p2 ,p3);

    MyCallback callback;

    bool const test = geometry::contacts_sphere_tetrahedron(sphere, tetrahedron, callback, false);

    BOOST_CHECK(!test);

    BOOST_CHECK_EQUAL(callback.m_contacts.size(), 0u);
  }
  {
    V const center = V::make( 1.0, -1.0, 1.0);
    T const radius = 1.0;

    V const p0 = V::make(0.0, 0.0, 0.0);
    V const p1 = V::make(1.0, 0.0, 0.0);
    V const p2 = V::make(0.0, 1.0, 0.0);
    V const p3 = V::make(0.0, 0.0, 1.0);

    geometry::Sphere<V>      const sphere      = geometry::make_sphere(center, radius);
    geometry::Tetrahedron<V> const tetrahedron = geometry::make_tetrahedron(p0, p1, p2 ,p3);

    MyCallback callback;

    bool const test = geometry::contacts_sphere_tetrahedron(sphere, tetrahedron, callback, false);

    BOOST_CHECK(!test);

    BOOST_CHECK_EQUAL(callback.m_contacts.size(), 0u);
  }
  {
    V const center = V::make(-1.0, 1.0, 1.0);
    T const radius = 1.0;

    V const p0 = V::make(0.0, 0.0, 0.0);
    V const p1 = V::make(1.0, 0.0, 0.0);
    V const p2 = V::make(0.0, 1.0, 0.0);
    V const p3 = V::make(0.0, 0.0, 1.0);

    geometry::Sphere<V>      const sphere      = geometry::make_sphere(center, radius);
    geometry::Tetrahedron<V> const tetrahedron = geometry::make_tetrahedron(p0, p1, p2 ,p3);

    MyCallback callback;

    bool const test = geometry::contacts_sphere_tetrahedron(sphere, tetrahedron, callback, false);

    BOOST_CHECK(!test);

    BOOST_CHECK_EQUAL(callback.m_contacts.size(), 0u);
  }
  // Separation in face regions
  {
    V const center = V::make( 0.2, 0.2, -2.1);
    T const radius = 1.0;

    V const p0 = V::make(0.0, 0.0, 0.0);
    V const p1 = V::make(1.0, 0.0, 0.0);
    V const p2 = V::make(0.0, 1.0, 0.0);
    V const p3 = V::make(0.0, 0.0, 1.0);

    geometry::Sphere<V>      const sphere      = geometry::make_sphere(center, radius);
    geometry::Tetrahedron<V> const tetrahedron = geometry::make_tetrahedron(p0, p1, p2 ,p3);

    MyCallback callback;

    bool const test = geometry::contacts_sphere_tetrahedron(sphere, tetrahedron, callback, false);

    BOOST_CHECK(!test);

    BOOST_CHECK_EQUAL(callback.m_contacts.size(), 0u);
  }
  {
    V const center = V::make( 0.2, -2.1,  0.2);
    T const radius = 1.0;

    V const p0 = V::make(0.0, 0.0, 0.0);
    V const p1 = V::make(1.0, 0.0, 0.0);
    V const p2 = V::make(0.0, 1.0, 0.0);
    V const p3 = V::make(0.0, 0.0, 1.0);

    geometry::Sphere<V>      const sphere      = geometry::make_sphere(center, radius);
    geometry::Tetrahedron<V> const tetrahedron = geometry::make_tetrahedron(p0, p1, p2 ,p3);

    MyCallback callback;

    bool const test = geometry::contacts_sphere_tetrahedron(sphere, tetrahedron, callback, false);

    BOOST_CHECK(!test);

    BOOST_CHECK_EQUAL(callback.m_contacts.size(), 0u);
  }
  {
    V const center = V::make( -2.1, 0.2,  0.2);
    T const radius = 1.0;

    V const p0 = V::make(0.0, 0.0, 0.0);
    V const p1 = V::make(1.0, 0.0, 0.0);
    V const p2 = V::make(0.0, 1.0, 0.0);
    V const p3 = V::make(0.0, 0.0, 1.0);

    geometry::Sphere<V>      const sphere      = geometry::make_sphere(center, radius);
    geometry::Tetrahedron<V> const tetrahedron = geometry::make_tetrahedron(p0, p1, p2 ,p3);

    MyCallback callback;

    bool const test = geometry::contacts_sphere_tetrahedron(sphere, tetrahedron, callback, false);

    BOOST_CHECK(!test);

    BOOST_CHECK_EQUAL(callback.m_contacts.size(), 0u);
  }
  {
    V const center = V::make( 2.0, 2.0, 2.0);
    T const radius = 1.0;

    V const p0 = V::make(0.0, 0.0, 0.0);
    V const p1 = V::make(1.0, 0.0, 0.0);
    V const p2 = V::make(0.0, 1.0, 0.0);
    V const p3 = V::make(0.0, 0.0, 1.0);

    geometry::Sphere<V>      const sphere      = geometry::make_sphere(center, radius);
    geometry::Tetrahedron<V> const tetrahedron = geometry::make_tetrahedron(p0, p1, p2 ,p3);

    MyCallback callback;

    bool const test = geometry::contacts_sphere_tetrahedron(sphere, tetrahedron, callback, false);

    BOOST_CHECK(!test);

    BOOST_CHECK_EQUAL(callback.m_contacts.size(), 0u);
  }
  // Contacts in vertex regions
  {
    V const center = V::make(1.9, 0.0, 0.0);
    T const radius = 1.0;

    V const p0 = V::make(0.0, 0.0, 0.0);
    V const p1 = V::make(1.0, 0.0, 0.0);
    V const p2 = V::make(0.0, 1.0, 0.0);
    V const p3 = V::make(0.0, 0.0, 1.0);

    geometry::Sphere<V>      const sphere      = geometry::make_sphere(center, radius);
    geometry::Tetrahedron<V> const tetrahedron = geometry::make_tetrahedron(p0, p1, p2 ,p3);

    MyCallback callback;

    bool const test = geometry::contacts_sphere_tetrahedron(sphere, tetrahedron, callback, false);

    BOOST_CHECK(test);

    BOOST_CHECK_EQUAL(callback.m_contacts.size(), 1u);

    BOOST_CHECK_CLOSE(callback.m_contacts[0].m_point(0), 1.0, 0.01);
    BOOST_CHECK_CLOSE(callback.m_contacts[0].m_point(1), 0.0, 0.01);
    BOOST_CHECK_CLOSE(callback.m_contacts[0].m_point(2), 0.0, 0.01);

    BOOST_CHECK_CLOSE(callback.m_contacts[0].m_normal(0),-1.0, 0.01);
    BOOST_CHECK_CLOSE(callback.m_contacts[0].m_normal(1), 0.0, 0.01);
    BOOST_CHECK_CLOSE(callback.m_contacts[0].m_normal(2), 0.0, 0.01);

    BOOST_CHECK_CLOSE(callback.m_contacts[0].m_distance, -0.1, 0.01);
  }
  {
    V const center = V::make(1.9, 0.0, 0.0);
    T const radius = 1.0;

    V const p0 = V::make(0.0, 0.0, 0.0);
    V const p1 = V::make(1.0, 0.0, 0.0);
    V const p2 = V::make(0.0, 1.0, 0.0);
    V const p3 = V::make(0.0, 0.0, 1.0);

    geometry::Sphere<V>      const sphere      = geometry::make_sphere(center, radius);
    geometry::Tetrahedron<V> const tetrahedron = geometry::make_tetrahedron(p0, p1, p2 ,p3);

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
    V const center = V::make(0.0, 1.9, 0.0);
    T const radius = 1.0;

    V const p0 = V::make(0.0, 0.0, 0.0);
    V const p1 = V::make(1.0, 0.0, 0.0);
    V const p2 = V::make(0.0, 1.0, 0.0);
    V const p3 = V::make(0.0, 0.0, 1.0);

    geometry::Sphere<V>      const sphere      = geometry::make_sphere(center, radius);
    geometry::Tetrahedron<V> const tetrahedron = geometry::make_tetrahedron(p0, p1, p2 ,p3);

    MyCallback callback;

    bool const test = geometry::contacts_sphere_tetrahedron(sphere, tetrahedron, callback, false);

    BOOST_CHECK(test);

    BOOST_CHECK_EQUAL(callback.m_contacts.size(), 1u);
  }
  {
    V const center = V::make(0.0, 0.0, 1.9);
    T const radius = 1.0;

    V const p0 = V::make(0.0, 0.0, 0.0);
    V const p1 = V::make(1.0, 0.0, 0.0);
    V const p2 = V::make(0.0, 1.0, 0.0);
    V const p3 = V::make(0.0, 0.0, 1.0);

    geometry::Sphere<V>      const sphere      = geometry::make_sphere(center, radius);
    geometry::Tetrahedron<V> const tetrahedron = geometry::make_tetrahedron(p0, p1, p2 ,p3);

    MyCallback callback;

    bool const test = geometry::contacts_sphere_tetrahedron(sphere, tetrahedron, callback, false);

    BOOST_CHECK(test);

    BOOST_CHECK_EQUAL(callback.m_contacts.size(), 1u);
  }
  {
    V const center = V::make(-0.5, -0.5, -0.5);
    T const radius = 1.0;

    V const p0 = V::make(0.0, 0.0, 0.0);
    V const p1 = V::make(1.0, 0.0, 0.0);
    V const p2 = V::make(0.0, 1.0, 0.0);
    V const p3 = V::make(0.0, 0.0, 1.0);

    geometry::Sphere<V>      const sphere      = geometry::make_sphere(center, radius);
    geometry::Tetrahedron<V> const tetrahedron = geometry::make_tetrahedron(p0, p1, p2 ,p3);

    MyCallback callback;

    bool const test = geometry::contacts_sphere_tetrahedron(sphere, tetrahedron, callback, false);

    BOOST_CHECK(test);

    BOOST_CHECK_EQUAL(callback.m_contacts.size(), 1u);
  }
  // Contact in edge regions
  {
    V const center = V::make(0.5, -0.5, -0.5);
    T const radius = 1.0;

    V const p0 = V::make(0.0, 0.0, 0.0);
    V const p1 = V::make(1.0, 0.0, 0.0);
    V const p2 = V::make(0.0, 1.0, 0.0);
    V const p3 = V::make(0.0, 0.0, 1.0);

    geometry::Sphere<V>      const sphere      = geometry::make_sphere(center, radius);
    geometry::Tetrahedron<V> const tetrahedron = geometry::make_tetrahedron(p0, p1, p2 ,p3);

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
    V const center = V::make(0.5, -0.5, -0.5);
    T const radius = 1.0;

    V const p0 = V::make(0.0, 0.0, 0.0);
    V const p1 = V::make(1.0, 0.0, 0.0);
    V const p2 = V::make(0.0, 1.0, 0.0);
    V const p3 = V::make(0.0, 0.0, 1.0);

    geometry::Sphere<V>      const sphere      = geometry::make_sphere(center, radius);
    geometry::Tetrahedron<V> const tetrahedron = geometry::make_tetrahedron(p0, p1, p2 ,p3);

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
    V const center = V::make(-0.5, 0.5, -0.5);
    T const radius = 1.0;

    V const p0 = V::make(0.0, 0.0, 0.0);
    V const p1 = V::make(1.0, 0.0, 0.0);
    V const p2 = V::make(0.0, 1.0, 0.0);
    V const p3 = V::make(0.0, 0.0, 1.0);

    geometry::Sphere<V>      const sphere      = geometry::make_sphere(center, radius);
    geometry::Tetrahedron<V> const tetrahedron = geometry::make_tetrahedron(p0, p1, p2 ,p3);

    MyCallback callback;

    bool const test = geometry::contacts_sphere_tetrahedron(sphere, tetrahedron, callback, false);

    BOOST_CHECK(test);

    BOOST_CHECK_EQUAL(callback.m_contacts.size(), 1u);
  }
  {
    V const center = V::make(-0.5, -0.5, 0.5);
    T const radius = 1.0;

    V const p0 = V::make(0.0, 0.0, 0.0);
    V const p1 = V::make(1.0, 0.0, 0.0);
    V const p2 = V::make(0.0, 1.0, 0.0);
    V const p3 = V::make(0.0, 0.0, 1.0);

    geometry::Sphere<V>      const sphere      = geometry::make_sphere(center, radius);
    geometry::Tetrahedron<V> const tetrahedron = geometry::make_tetrahedron(p0, p1, p2 ,p3);

    MyCallback callback;

    bool const test = geometry::contacts_sphere_tetrahedron(sphere, tetrahedron,  callback, false);

    BOOST_CHECK(test);

    BOOST_CHECK_EQUAL(callback.m_contacts.size(), 1u);
  }
  {
    V const center = V::make( 0.5, 0.5, -0.5);
    T const radius = 1.0;

    V const p0 = V::make(0.0, 0.0, 0.0);
    V const p1 = V::make(1.0, 0.0, 0.0);
    V const p2 = V::make(0.0, 1.0, 0.0);
    V const p3 = V::make(0.0, 0.0, 1.0);

    geometry::Sphere<V>      const sphere      = geometry::make_sphere(center, radius);
    geometry::Tetrahedron<V> const tetrahedron = geometry::make_tetrahedron(p0, p1, p2 ,p3);

    MyCallback callback;

    bool const test = geometry::contacts_sphere_tetrahedron(sphere, tetrahedron, callback, false);

    BOOST_CHECK(test);

    BOOST_CHECK_EQUAL(callback.m_contacts.size(), 1u);
  }
  {
    V const center = V::make( 0.5, -0.5, 0.5);
    T const radius = 1.0;

    V const p0 = V::make(0.0, 0.0, 0.0);
    V const p1 = V::make(1.0, 0.0, 0.0);
    V const p2 = V::make(0.0, 1.0, 0.0);
    V const p3 = V::make(0.0, 0.0, 1.0);

    geometry::Sphere<V>      const sphere      = geometry::make_sphere(center, radius);
    geometry::Tetrahedron<V> const tetrahedron = geometry::make_tetrahedron(p0, p1, p2 ,p3);

    MyCallback callback;

    bool const test = geometry::contacts_sphere_tetrahedron(sphere, tetrahedron, callback, false);

    BOOST_CHECK(test);

    BOOST_CHECK_EQUAL(callback.m_contacts.size(), 1u);
  }
  {
    V const center = V::make(-0.5, 0.5, 0.5);
    T const radius = 1.0;

    V const p0 = V::make(0.0, 0.0, 0.0);
    V const p1 = V::make(1.0, 0.0, 0.0);
    V const p2 = V::make(0.0, 1.0, 0.0);
    V const p3 = V::make(0.0, 0.0, 1.0);

    geometry::Sphere<V>      const sphere      = geometry::make_sphere(center, radius);
    geometry::Tetrahedron<V> const tetrahedron = geometry::make_tetrahedron(p0, p1, p2 ,p3);

    MyCallback callback;

    bool const test = geometry::contacts_sphere_tetrahedron(sphere, tetrahedron, callback, false);

    BOOST_CHECK(test);

    BOOST_CHECK_EQUAL(callback.m_contacts.size(), 1u);
  }
  // Contact in face regions
  {
    V const center = V::make( 0.2, 0.2, -0.9);
    T const radius = 1.0;

    V const p0 = V::make(0.0, 0.0, 0.0);
    V const p1 = V::make(1.0, 0.0, 0.0);
    V const p2 = V::make(0.0, 1.0, 0.0);
    V const p3 = V::make(0.0, 0.0, 1.0);

    geometry::Sphere<V>      const sphere      = geometry::make_sphere(center, radius);
    geometry::Tetrahedron<V> const tetrahedron = geometry::make_tetrahedron(p0, p1, p2 ,p3);

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
    V const center = V::make( 0.2, 0.2, -0.9);
    T const radius = 1.0;

    V const p0 = V::make(0.0, 0.0, 0.0);
    V const p1 = V::make(1.0, 0.0, 0.0);
    V const p2 = V::make(0.0, 1.0, 0.0);
    V const p3 = V::make(0.0, 0.0, 1.0);

    geometry::Sphere<V>      const sphere      = geometry::make_sphere(center, radius);
    geometry::Tetrahedron<V> const tetrahedron = geometry::make_tetrahedron(p0, p1, p2 ,p3);

    MyCallback callback;

    bool const test = geometry::contacts_sphere_tetrahedron(sphere, tetrahedron, callback, true);

    BOOST_CHECK(test);

    BOOST_CHECK_EQUAL(callback.m_contacts.size(), 1u);

    BOOST_CHECK_CLOSE(callback.m_contacts[0].m_point(0), 0.2, 0.01);
    BOOST_CHECK_CLOSE(callback.m_contacts[0].m_point(1), 0.2, 0.01);
    BOOST_CHECK_CLOSE(callback.m_contacts[0].m_point(2), 0.0, 0.01);

    BOOST_CHECK_CLOSE(callback.m_contacts[0].m_normal(0), 0.0, 0.01);
    BOOST_CHECK_CLOSE(callback.m_contacts[0].m_normal(1), 0.0, 0.01);
    BOOST_CHECK_CLOSE(callback.m_contacts[0].m_normal(2),-1.0, 0.01);

    BOOST_CHECK_CLOSE(callback.m_contacts[0].m_distance, -0.1, 0.01);
  }
  {
    V const center = V::make( 0.2, -0.9,  0.2);
    T const radius = 1.0;

    V const p0 = V::make(0.0, 0.0, 0.0);
    V const p1 = V::make(1.0, 0.0, 0.0);
    V const p2 = V::make(0.0, 1.0, 0.0);
    V const p3 = V::make(0.0, 0.0, 1.0);

    geometry::Sphere<V>      const sphere      = geometry::make_sphere(center, radius);
    geometry::Tetrahedron<V> const tetrahedron = geometry::make_tetrahedron(p0, p1, p2 ,p3);

    MyCallback callback;

    bool const test = geometry::contacts_sphere_tetrahedron(sphere, tetrahedron, callback, false);

    BOOST_CHECK(test);

    BOOST_CHECK_EQUAL(callback.m_contacts.size(), 1u);
  }
  {
    V const center = V::make( -0.9, 0.2,  0.2);
    T const radius = 1.0;

    V const p0 = V::make(0.0, 0.0, 0.0);
    V const p1 = V::make(1.0, 0.0, 0.0);
    V const p2 = V::make(0.0, 1.0, 0.0);
    V const p3 = V::make(0.0, 0.0, 1.0);

    geometry::Sphere<V>      const sphere      = geometry::make_sphere(center, radius);
    geometry::Tetrahedron<V> const tetrahedron = geometry::make_tetrahedron(p0, p1, p2 ,p3);

    MyCallback callback;

    bool const test = geometry::contacts_sphere_tetrahedron(sphere, tetrahedron, callback, false);

    BOOST_CHECK(test);

    BOOST_CHECK_EQUAL(callback.m_contacts.size(), 1u);
  }
  {
    V const center = V::make( 0.7, 0.7, 0.7);
    T const radius = 1.0;

    V const p0 = V::make(0.0, 0.0, 0.0);
    V const p1 = V::make(1.0, 0.0, 0.0);
    V const p2 = V::make(0.0, 1.0, 0.0);
    V const p3 = V::make(0.0, 0.0, 1.0);

    geometry::Sphere<V>      const sphere      = geometry::make_sphere(center, radius);
    geometry::Tetrahedron<V> const tetrahedron = geometry::make_tetrahedron(p0, p1, p2 ,p3);
    
    MyCallback callback;
    
    bool const test = geometry::contacts_sphere_tetrahedron(sphere, tetrahedron, callback, false);
    
    BOOST_CHECK(test);

    BOOST_CHECK_EQUAL(callback.m_contacts.size(), 1u);
  }

  // sphere inside
  {
    V const center = V::make( 0.3, 0.2, 0.1);
    T const radius = 1.0;

    V const p0 = V::make(0.0, 0.0, 0.0);
    V const p1 = V::make(1.0, 0.0, 0.0);
    V const p2 = V::make(0.0, 1.0, 0.0);
    V const p3 = V::make(0.0, 0.0, 1.0);

    geometry::Sphere<V>      const sphere      = geometry::make_sphere(center, radius);
    geometry::Tetrahedron<V> const tetrahedron = geometry::make_tetrahedron(p0, p1, p2 ,p3);

    MyCallback callback;

    bool const test = geometry::contacts_sphere_tetrahedron(sphere, tetrahedron, callback, false);

    BOOST_CHECK(test);

    BOOST_CHECK_EQUAL(callback.m_contacts.size(), 1u);
  }

}

BOOST_AUTO_TEST_SUITE_END();
