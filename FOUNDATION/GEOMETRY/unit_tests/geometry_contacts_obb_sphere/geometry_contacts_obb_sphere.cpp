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

BOOST_AUTO_TEST_CASE(contacts_obb_sphere_test)
{
  // Touching rigth side
  {
    V const centerA   = V::make(0.0, 0.0, 0.0);
    V const half_extA = V::make(1.0, 1.0, 1.0);
    Q const qA        = Q::identity();

    V const centerB = V::make(2.0, 0.0, 0.0);
    T const radiusB = 1.0;

    geometry::OBB<MT>   const obb    = geometry::make_obb<MT>(centerA, qA, half_extA);
    geometry::Sphere<V> const sphere = geometry::make_sphere(centerB, radiusB);

    MyCallback callback;

    geometry::contacts_obb_sphere(obb, sphere, 0.0, callback, false);

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
    V const centerA   = V::make(0.0, 0.0, 0.0);
    V const half_extA = V::make(1.0, 1.0, 1.0);
    Q const qA        = Q::identity();

    V const centerB = V::make(2.1, 0.0, 0.0);
    T const radiusB = 1.0;

    geometry::OBB<MT>   const obb    = geometry::make_obb<MT>(centerA, qA, half_extA);
    geometry::Sphere<V> const sphere = geometry::make_sphere(centerB, radiusB);

    MyCallback callback;

    geometry::contacts_obb_sphere(obb, sphere, 0.0, callback, false);

    BOOST_CHECK(callback.m_contacts.size() == 0u);
  }
  // Penetration rigth side
  {
    V const centerA   = V::make(0.0, 0.0, 0.0);
    V const half_extA = V::make(1.0, 1.0, 1.0);
    Q const qA        = Q::identity();

    V const centerB = V::make(1.5, 0.0, 0.0);
    T const radiusB = 1.0;

    geometry::OBB<MT>   const obb    = geometry::make_obb<MT>(centerA, qA, half_extA);
    geometry::Sphere<V> const sphere = geometry::make_sphere(centerB, radiusB);

    MyCallback callback;

    geometry::contacts_obb_sphere(obb, sphere, 0.0, callback, false);

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
    V const centerA   = V::make(0.0, 0.0, 0.0);
    V const half_extA = V::make(1.0, 1.0, 1.0);
    Q const qA        = Q::identity();

    V const centerB = V::make(1.5, 0.0, 0.0);
    T const radiusB = 1.0;

    geometry::OBB<MT>   const obb    = geometry::make_obb<MT>(centerA, qA, half_extA);
    geometry::Sphere<V> const sphere = geometry::make_sphere(centerB, radiusB);

    MyCallback callback;

    geometry::contacts_obb_sphere(obb, sphere, 0.0, callback, true);

    BOOST_CHECK(callback.m_contacts.size() == 1u);

    BOOST_CHECK_CLOSE(callback.m_contacts[0].m_point[0], 1.0, 0.01);
    BOOST_CHECK_CLOSE(callback.m_contacts[0].m_point[1], 0.0, 0.01);
    BOOST_CHECK_CLOSE(callback.m_contacts[0].m_point[2], 0.0, 0.01);

    BOOST_CHECK_CLOSE(callback.m_contacts[0].m_normal[0],-1.0, 0.01);
    BOOST_CHECK_CLOSE(callback.m_contacts[0].m_normal[1], 0.0, 0.01);
    BOOST_CHECK_CLOSE(callback.m_contacts[0].m_normal[2], 0.0, 0.01);

    BOOST_CHECK_CLOSE(callback.m_contacts[0].m_distance, -0.5, 0.01);
    
  }

  // Touching left side
  {
    V const centerA   = V::make(0.0, 0.0, 0.0);
    V const half_extA = V::make(1.0, 1.0, 1.0);
    Q const qA        = Q::identity();

    V const centerB = V::make(-2.0, 0.0, 0.0);
    T const radiusB = 1.0;

    geometry::OBB<MT>   const obb    = geometry::make_obb<MT>(centerA, qA, half_extA);
    geometry::Sphere<V> const sphere = geometry::make_sphere(centerB, radiusB);

    MyCallback callback;

    geometry::contacts_obb_sphere(obb, sphere, 0.0, callback, false);

    BOOST_CHECK(callback.m_contacts.size() == 1u);

    BOOST_CHECK_CLOSE(callback.m_contacts[0].m_point[0],-1.0, 0.01);
    BOOST_CHECK_CLOSE(callback.m_contacts[0].m_point[1], 0.0, 0.01);
    BOOST_CHECK_CLOSE(callback.m_contacts[0].m_point[2], 0.0, 0.01);

    BOOST_CHECK_CLOSE(callback.m_contacts[0].m_normal[0],-1.0, 0.01);
    BOOST_CHECK_CLOSE(callback.m_contacts[0].m_normal[1], 0.0, 0.01);
    BOOST_CHECK_CLOSE(callback.m_contacts[0].m_normal[2], 0.0, 0.01);

    BOOST_CHECK_CLOSE(callback.m_contacts[0].m_distance, 0.0, 0.01);

  }
  // Separating left side
  {
    V const centerA   = V::make(0.0, 0.0, 0.0);
    V const half_extA = V::make(1.0, 1.0, 1.0);
    Q const qA        = Q::identity();

    V const centerB = V::make(-2.1, 0.0, 0.0);
    T const radiusB = 1.0;

    geometry::OBB<MT>   const obb    = geometry::make_obb<MT>(centerA, qA, half_extA);
    geometry::Sphere<V> const sphere = geometry::make_sphere(centerB, radiusB);

    MyCallback callback;

    geometry::contacts_obb_sphere(obb, sphere, 0.0, callback, false);

    BOOST_CHECK(callback.m_contacts.size() == 0u);
  }
  // Penetration left side
  {
    V const centerA   = V::make(0.0, 0.0, 0.0);
    V const half_extA = V::make(1.0, 1.0, 1.0);
    Q const qA        = Q::identity();

    V const centerB = V::make(-1.5, 0.0, 0.0);
    T const radiusB = 1.0;

    geometry::OBB<MT>   const obb    = geometry::make_obb<MT>(centerA, qA, half_extA);
    geometry::Sphere<V> const sphere = geometry::make_sphere(centerB, radiusB);

    MyCallback callback;

    geometry::contacts_obb_sphere(obb, sphere, 0.0, callback, false);

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
    V const centerA   = V::make(0.0, 0.0, 0.0);
    V const half_extA = V::make(1.0, 1.0, 1.0);
    Q const qA        = Q::identity();

    V const centerB = V::make(-1.5, 0.0, 0.0);
    T const radiusB = 1.0;

    geometry::OBB<MT>   const obb    = geometry::make_obb<MT>(centerA, qA, half_extA);
    geometry::Sphere<V> const sphere = geometry::make_sphere(centerB, radiusB);

    MyCallback callback;

    geometry::contacts_obb_sphere(obb, sphere, 0.0, callback, true);

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
