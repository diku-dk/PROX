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

BOOST_AUTO_TEST_CASE(contacts_sphere_sphere_test)
{

  // Penetration
  {
    V const centerA = V::make(0.0, 0.0, 0.0);
    T const radiusA = 2.0;
    V const centerB = V::make(5.0, 0.0, 0.0);
    T const radiusB = 4.0;

    geometry::Sphere<V> const A = geometry::make_sphere(centerA, radiusA);
    geometry::Sphere<V> const B = geometry::make_sphere(centerB, radiusB);

    MyCallback callback;

    geometry::contacts_sphere_sphere(A, B, 0.0, callback);

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
    V const centerA = V::make(0.0, 0.0, 0.0);
    T const radiusA = 2.0;
    V const centerB = V::make(6.0, 0.0, 0.0);
    T const radiusB = 4.0;

    geometry::Sphere<V> const A = geometry::make_sphere(centerA, radiusA);
    geometry::Sphere<V> const B = geometry::make_sphere(centerB, radiusB);

    MyCallback callback;

    geometry::contacts_sphere_sphere(A, B, 0.0, callback);

    BOOST_CHECK(callback.m_contacts.size() == 1u);

    BOOST_CHECK_CLOSE(callback.m_contacts[0].m_point[0], 2.0, 0.01);
    BOOST_CHECK_CLOSE(callback.m_contacts[0].m_point[1], 0.0, 0.01);
    BOOST_CHECK_CLOSE(callback.m_contacts[0].m_point[2], 0.0, 0.01);

    BOOST_CHECK_CLOSE(callback.m_contacts[0].m_normal[0], 1.0, 0.01);
    BOOST_CHECK_CLOSE(callback.m_contacts[0].m_normal[1], 0.0, 0.01);
    BOOST_CHECK_CLOSE(callback.m_contacts[0].m_normal[2], 0.0, 0.01);

    BOOST_CHECK_CLOSE(callback.m_contacts[0].m_distance,  0.0, 0.01);
    
    
  }

  // Separation
  {
    V const centerA = V::make(0.0, 0.0, 0.0);
    T const radiusA = 2.0;
    V const centerB = V::make(7.0, 0.0, 0.0);
    T const radiusB = 4.0;

    geometry::Sphere<V> const A = geometry::make_sphere(centerA, radiusA);
    geometry::Sphere<V> const B = geometry::make_sphere(centerB, radiusB);

    MyCallback callback;

    geometry::contacts_sphere_sphere(A, B, 0.0, callback);

    BOOST_CHECK(callback.m_contacts.size() == 0u);
  }

}

BOOST_AUTO_TEST_SUITE_END();
