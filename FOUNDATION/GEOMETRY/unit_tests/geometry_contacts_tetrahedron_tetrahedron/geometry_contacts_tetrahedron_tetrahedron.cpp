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

BOOST_AUTO_TEST_CASE(contacts_tetrahedron_tetrahedron_test)
{
  // B tip inside left-side of A
  {
    V const a0 = V::make(0.0, 0.0, 0.0);
    V const a1 = V::make(1.0, 0.0, 0.0);
    V const a2 = V::make(0.0, 1.0, 0.0);
    V const a3 = V::make(0.0, 0.0, 1.0);

    V const offset = V::make(-0.7, 0.2, 0.2);

    V const b0 = V::make(0.0, 0.0, 0.0) + offset;
    V const b1 = V::make(1.0, 0.0, 0.0) + offset;
    V const b2 = V::make(0.0, 1.0, 0.0) + offset;
    V const b3 = V::make(0.0, 0.0, 1.0) + offset;

    geometry::Tetrahedron<V> const A = geometry::make_tetrahedron(a0, a1, a2 ,a3);
    geometry::Tetrahedron<V> const B = geometry::make_tetrahedron(b0, b1, b2 ,b3);

    MyCallback callback;

    bool const test = geometry::contacts_tetrahedron_tetrahedron(A, B, callback );

    BOOST_CHECK(test);

    BOOST_CHECK_EQUAL(callback.m_contacts.size(), 3u);

    BOOST_CHECK_CLOSE( callback.m_contacts[0].m_point(0),   0.15, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[0].m_point(1),   0.2,  0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[0].m_point(2),   0.2,  0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[0].m_normal(0), -1.0,  0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[0].m_normal(1),  0.0,  0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[0].m_normal(2),  0.0,  0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[0].m_distance,  -0.3,  0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[1].m_point(0),   0.15, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[1].m_point(1),   0.5,  0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[1].m_point(2),   0.2,  0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[1].m_normal(0), -1.0,  0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[1].m_normal(1),  0.0,  0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[1].m_normal(2),  0.0,  0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[1].m_distance,  -0.3,  0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[2].m_point(0),   0.15, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[2].m_point(1),   0.2,  0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[2].m_point(2),   0.5,  0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[2].m_normal(0), -1.0,  0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[2].m_normal(1),  0.0,  0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[2].m_normal(2),  0.0,  0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[2].m_distance,  -0.3,  0.01);

  }
  // B tip inside front-side of A
  {
    V const a0 = V::make(0.0, 0.0, 0.0);
    V const a1 = V::make(1.0, 0.0, 0.0);
    V const a2 = V::make(0.0, 1.0, 0.0);
    V const a3 = V::make(0.0, 0.0, 1.0);

    V const offset = V::make(0.3, -0.7, 0.2);

    V const b0 = V::make(0.0, 0.0, 0.0) + offset;
    V const b1 = V::make(1.0, 0.0, 0.0) + offset;
    V const b2 = V::make(0.0, 1.0, 0.0) + offset;
    V const b3 = V::make(0.0, 0.0, 1.0) + offset;

    geometry::Tetrahedron<V> const A = geometry::make_tetrahedron(a0, a1, a2 ,a3);
    geometry::Tetrahedron<V> const B = geometry::make_tetrahedron(b0, b1, b2 ,b3);

    MyCallback callback;

    bool const test = geometry::contacts_tetrahedron_tetrahedron(A, B, callback );

    BOOST_CHECK(test);

    BOOST_CHECK_EQUAL(callback.m_contacts.size(), 3u);

    BOOST_CHECK_CLOSE( callback.m_contacts[0].m_point(0),   0.30000001192, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[0].m_point(1),   0.15000000596, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[0].m_point(2),   0.20000000298, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[0].m_normal(0),  0.0,           0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[0].m_normal(1), -1.0,           0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[0].m_normal(2),  0.0,           0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[0].m_distance,  -0.30000001192, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[1].m_point(0),   0.60000002384, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[1].m_point(1),   0.15000000596, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[1].m_point(2),   0.20000000298, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[1].m_normal(0),  0.0,           0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[1].m_normal(1), -1.0,           0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[1].m_normal(2),  0.0,           0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[1].m_distance,  -0.30000001192, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[2].m_point(0),   0.30000001192, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[2].m_point(1),   0.15000000596, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[2].m_point(2),   0.5,           0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[2].m_normal(0),  0.0,           0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[2].m_normal(1), -1.0,           0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[2].m_normal(2),  0.0,           0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[2].m_distance,  -0.30000001192, 0.01);
  }
  // B tip inside bottom-side of A
  {
    V const a0 = V::make(0.0, 0.0, 0.0);
    V const a1 = V::make(1.0, 0.0, 0.0);
    V const a2 = V::make(0.0, 1.0, 0.0);
    V const a3 = V::make(0.0, 0.0, 1.0);

    V const offset = V::make(0.3, 0.3, -0.7);

    V const b0 = V::make(0.0, 0.0, 0.0) + offset;
    V const b1 = V::make(1.0, 0.0, 0.0) + offset;
    V const b2 = V::make(0.0, 1.0, 0.0) + offset;
    V const b3 = V::make(0.0, 0.0, 1.0) + offset;

    geometry::Tetrahedron<V> const A = geometry::make_tetrahedron(a0, a1, a2 ,a3);
    geometry::Tetrahedron<V> const B = geometry::make_tetrahedron(b0, b1, b2 ,b3);

    MyCallback callback;

    bool const test = geometry::contacts_tetrahedron_tetrahedron(A, B, callback );

    BOOST_CHECK(test);

    BOOST_CHECK_EQUAL(callback.m_contacts.size(), 4u);

    BOOST_CHECK_CLOSE( callback.m_contacts[0].m_point(0),   0.3749999702,  0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[0].m_point(1),   0.3749999702,  0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[0].m_point(2),   0.30000001192, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[0].m_normal(0),  0.70710676908, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[0].m_normal(1),  0.70710676908, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[0].m_normal(2),  0.0,           0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[0].m_distance,  -0.21213203669, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[1].m_point(0),   0.3749999702,  0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[1].m_point(1),   0.3749999702,  0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[1].m_point(2),   0.0,           0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[1].m_normal(0),  0.70710676908, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[1].m_normal(1),  0.70710676908, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[1].m_normal(2),  0.0,           0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[1].m_distance,  -0.21213203669, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[2].m_point(0),   0.52499997616, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[2].m_point(1),   0.22499996424, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[2].m_point(2),   0.0,           0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[2].m_normal(0),  0.70710676908, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[2].m_normal(1),  0.70710676908, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[2].m_normal(2),  0.0,           0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[2].m_distance,  -0.21213203669, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[3].m_point(0),   0.22499996424, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[3].m_point(1),   0.52499997616, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[3].m_point(2),   0.0,           0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[3].m_normal(0),  0.70710676908, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[3].m_normal(1),  0.70710676908, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[3].m_normal(2),  0.0,           0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[3].m_distance,  -0.21213203669, 0.01);

  }
  // B tip inside oblique-side of A
  {
    V const a0 = V::make(0.0, 0.0, 0.0);
    V const a1 = V::make(1.0, 0.0, 0.0);
    V const a2 = V::make(0.0, 1.0, 0.0);
    V const a3 = V::make(0.0, 0.0, 1.0);

    V const offset = V::make(0.3, 0.3, 0.3);

    V const b0 = V::make(0.0, 0.0, 0.0) + offset;
    V const b1 = V::make(1.0, 0.0, 0.0) + offset;
    V const b2 = V::make(0.0, 1.0, 0.0) + offset;
    V const b3 = V::make(0.0, 0.0, 1.0) + offset;

    geometry::Tetrahedron<V> const A = geometry::make_tetrahedron(a0, a1, a2 ,a3);
    geometry::Tetrahedron<V> const B = geometry::make_tetrahedron(b0, b1, b2 ,b3);

    MyCallback callback;

    bool const test = geometry::contacts_tetrahedron_tetrahedron(A, B, callback );

    BOOST_CHECK(test);

    BOOST_CHECK_EQUAL(callback.m_contacts.size(), 4u);

    BOOST_CHECK_CLOSE( callback.m_contacts[0].m_point(0), 0.316667, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[0].m_point(1), 0.316667, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[0].m_point(2), 0.316667, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[0].m_normal(0), 0.57735, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[0].m_normal(1), 0.57735, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[0].m_normal(2), 0.57735, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[0].m_distance, -0.057735, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[1].m_point(0), 0.383333, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[1].m_point(1), 0.283333, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[1].m_point(2), 0.283333, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[1].m_normal(0), 0.57735, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[1].m_normal(1), 0.57735, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[1].m_normal(2), 0.57735, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[1].m_distance, -0.057735, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[2].m_point(0), 0.283333, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[2].m_point(1), 0.383333, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[2].m_point(2), 0.283333, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[2].m_normal(0), 0.57735, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[2].m_normal(1), 0.57735, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[2].m_normal(2), 0.57735, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[2].m_distance, -0.057735, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[3].m_point(0), 0.283333, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[3].m_point(1), 0.283333, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[3].m_point(2), 0.383333, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[3].m_normal(0), 0.57735, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[3].m_normal(1), 0.57735, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[3].m_normal(2), 0.57735, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[3].m_distance, -0.057735, 0.01);
  }
  // B inside of A
  {
    V const a0 = V::make(0.0, 0.0, 0.0);
    V const a1 = V::make(1.0, 0.0, 0.0);
    V const a2 = V::make(0.0, 1.0, 0.0);
    V const a3 = V::make(0.0, 0.0, 1.0);

    V const b0 = V::make(0.1, 0.1, 0.1);
    V const b1 = V::make(0.9, 0.0, 0.0);
    V const b2 = V::make(0.0, 0.9, 0.0);
    V const b3 = V::make(0.0, 0.0, 0.9);

    geometry::Tetrahedron<V> const A = geometry::make_tetrahedron(a0, a1, a2 ,a3);
    geometry::Tetrahedron<V> const B = geometry::make_tetrahedron(b0, b1, b2 ,b3);

    MyCallback callback1;
    MyCallback callback2;

    bool const test1 = geometry::contacts_tetrahedron_tetrahedron(A, B, callback1 );
    bool const test2 = geometry::contacts_tetrahedron_tetrahedron(B, A, callback2 );

    BOOST_CHECK(test1);
    BOOST_CHECK(test2);

    BOOST_CHECK_EQUAL(callback1.m_contacts.size(), 4u);
    BOOST_CHECK_EQUAL(callback2.m_contacts.size(), 4u);


    BOOST_CHECK_CLOSE( callback1.m_contacts[0].m_point(0), 0.2, 0.01);
    BOOST_CHECK_CLOSE( callback1.m_contacts[0].m_point(1), 0.2, 0.01);
    BOOST_CHECK_CLOSE( callback1.m_contacts[0].m_point(2), 0.2, 0.01);
    BOOST_CHECK_CLOSE( callback1.m_contacts[0].m_normal(0), 0.57735, 0.01);
    BOOST_CHECK_CLOSE( callback1.m_contacts[0].m_normal(1), 0.57735, 0.01);
    BOOST_CHECK_CLOSE( callback1.m_contacts[0].m_normal(2), 0.57735, 0.01);
    BOOST_CHECK_CLOSE( callback1.m_contacts[0].m_distance, -0.34641, 0.01);
    BOOST_CHECK_CLOSE( callback1.m_contacts[1].m_point(0), 0.8, 0.01);
    BOOST_CHECK_CLOSE( callback1.m_contacts[1].m_point(1), -0.1, 0.01);
    BOOST_CHECK_CLOSE( callback1.m_contacts[1].m_point(2), -0.1, 0.01);
    BOOST_CHECK_CLOSE( callback1.m_contacts[1].m_normal(0), 0.57735, 0.01);
    BOOST_CHECK_CLOSE( callback1.m_contacts[1].m_normal(1), 0.57735, 0.01);
    BOOST_CHECK_CLOSE( callback1.m_contacts[1].m_normal(2), 0.57735, 0.01);
    BOOST_CHECK_CLOSE( callback1.m_contacts[1].m_distance, -0.34641, 0.01);
    BOOST_CHECK_CLOSE( callback1.m_contacts[2].m_point(0), -0.1, 0.01);
    BOOST_CHECK_CLOSE( callback1.m_contacts[2].m_point(1), 0.8, 0.01);
    BOOST_CHECK_CLOSE( callback1.m_contacts[2].m_point(2), -0.1, 0.01);
    BOOST_CHECK_CLOSE( callback1.m_contacts[2].m_normal(0), 0.57735, 0.01);
    BOOST_CHECK_CLOSE( callback1.m_contacts[2].m_normal(1), 0.57735, 0.01);
    BOOST_CHECK_CLOSE( callback1.m_contacts[2].m_normal(2), 0.57735, 0.01);
    BOOST_CHECK_CLOSE( callback1.m_contacts[2].m_distance, -0.34641, 0.01);
    BOOST_CHECK_CLOSE( callback1.m_contacts[3].m_point(0), -0.1, 0.01);
    BOOST_CHECK_CLOSE( callback1.m_contacts[3].m_point(1), -0.1, 0.01);
    BOOST_CHECK_CLOSE( callback1.m_contacts[3].m_point(2), 0.8, 0.01);
    BOOST_CHECK_CLOSE( callback1.m_contacts[3].m_normal(0), 0.57735, 0.01);
    BOOST_CHECK_CLOSE( callback1.m_contacts[3].m_normal(1), 0.57735, 0.01);
    BOOST_CHECK_CLOSE( callback1.m_contacts[3].m_normal(2), 0.57735, 0.01);
    BOOST_CHECK_CLOSE( callback1.m_contacts[3].m_distance, -0.34641, 0.01);

    BOOST_CHECK_CLOSE( callback2.m_contacts[0].m_point(0), 0.2, 0.01);
    BOOST_CHECK_CLOSE( callback2.m_contacts[0].m_point(1), 0.2, 0.01);
    BOOST_CHECK_CLOSE( callback2.m_contacts[0].m_point(2), 0.2, 0.01);
    BOOST_CHECK_CLOSE( callback2.m_contacts[0].m_normal(0), -0.57735, 0.01);
    BOOST_CHECK_CLOSE( callback2.m_contacts[0].m_normal(1), -0.57735, 0.01);
    BOOST_CHECK_CLOSE( callback2.m_contacts[0].m_normal(2), -0.57735, 0.01);
    BOOST_CHECK_CLOSE( callback2.m_contacts[0].m_distance, -0.34641, 0.01);
    BOOST_CHECK_CLOSE( callback2.m_contacts[1].m_point(0), 0.8, 0.01);
    BOOST_CHECK_CLOSE( callback2.m_contacts[1].m_point(1), -0.1, 0.01);
    BOOST_CHECK_CLOSE( callback2.m_contacts[1].m_point(2), -0.1, 0.01);
    BOOST_CHECK_CLOSE( callback2.m_contacts[1].m_normal(0), -0.57735, 0.01);
    BOOST_CHECK_CLOSE( callback2.m_contacts[1].m_normal(1), -0.57735, 0.01);
    BOOST_CHECK_CLOSE( callback2.m_contacts[1].m_normal(2), -0.57735, 0.01);
    BOOST_CHECK_CLOSE( callback2.m_contacts[1].m_distance, -0.34641, 0.01);
    BOOST_CHECK_CLOSE( callback2.m_contacts[2].m_point(0), -0.1, 0.01);
    BOOST_CHECK_CLOSE( callback2.m_contacts[2].m_point(1), 0.8, 0.01);
    BOOST_CHECK_CLOSE( callback2.m_contacts[2].m_point(2), -0.1, 0.01);
    BOOST_CHECK_CLOSE( callback2.m_contacts[2].m_normal(0), -0.57735, 0.01);
    BOOST_CHECK_CLOSE( callback2.m_contacts[2].m_normal(1), -0.57735, 0.01);
    BOOST_CHECK_CLOSE( callback2.m_contacts[2].m_normal(2), -0.57735, 0.01);
    BOOST_CHECK_CLOSE( callback2.m_contacts[2].m_distance, -0.34641, 0.01);
    BOOST_CHECK_CLOSE( callback2.m_contacts[3].m_point(0), -0.1, 0.01);
    BOOST_CHECK_CLOSE( callback2.m_contacts[3].m_point(1), -0.1, 0.01);
    BOOST_CHECK_CLOSE( callback2.m_contacts[3].m_point(2), 0.8, 0.01);
    BOOST_CHECK_CLOSE( callback2.m_contacts[3].m_normal(0), -0.57735, 0.01);
    BOOST_CHECK_CLOSE( callback2.m_contacts[3].m_normal(1), -0.57735, 0.01);
    BOOST_CHECK_CLOSE( callback2.m_contacts[3].m_normal(2), -0.57735, 0.01);
    BOOST_CHECK_CLOSE( callback2.m_contacts[3].m_distance, -0.34641, 0.01);
  }
  // Edge from B goes through A and vice-versa
  {
    V const a0 = V::make(0.0, 0.0, 0.0);
    V const a1 = V::make(1.0, 0.0, 0.0);
    V const a2 = V::make(0.0, 1.0, 0.0);
    V const a3 = V::make(0.0, 0.0, 1.0);

    V const b0 = V::make(0.3, 0.3, -1.0);
    V const b1 = V::make(1.3, 0.3,  0.0);
    V const b2 = V::make(0.3, 1.3,  0.0);
    V const b3 = V::make(0.3, 0.3,  1.0);

    geometry::Tetrahedron<V> const A = geometry::make_tetrahedron(a0, a1, a2 ,a3);
    geometry::Tetrahedron<V> const B = geometry::make_tetrahedron(b0, b1, b2 ,b3);

    MyCallback callback;

    bool const test = geometry::contacts_tetrahedron_tetrahedron(A, B, callback );

    BOOST_CHECK(test);

    BOOST_CHECK_EQUAL(callback.m_contacts.size(), 4u);

    BOOST_CHECK_CLOSE( callback.m_contacts[0].m_point(0), 0.4, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[0].m_point(1), 0.4, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[0].m_point(2), 0.4, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[0].m_normal(0), 0.707107, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[0].m_normal(1), 0.707107, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[0].m_normal(2), -0, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[0].m_distance, -0.282843, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[1].m_point(0), 0.4, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[1].m_point(1), 0.4, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[1].m_point(2), 0, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[1].m_normal(0), 0.707107, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[1].m_normal(1), 0.707107, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[1].m_normal(2), -0, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[1].m_distance, -0.282843, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[2].m_point(0), 0.2, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[2].m_point(1), 0.6, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[2].m_point(2), 0, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[2].m_normal(0), 0.707107, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[2].m_normal(1), 0.707107, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[2].m_normal(2), -0, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[2].m_distance, -0.282843, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[3].m_point(0), 0.6, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[3].m_point(1), 0.2, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[3].m_point(2), 0, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[3].m_normal(0), 0.707107, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[3].m_normal(1), 0.707107, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[3].m_normal(2), -0, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[3].m_distance, -0.282843, 0.01);

  }
  // A is impaled by B
  {
    V const a0 = V::make(0.0, 0.0, 0.0);
    V const a1 = V::make(1.0, 0.0, 0.0);
    V const a2 = V::make(0.0, 1.0, 0.0);
    V const a3 = V::make(0.0, 0.0, 1.0);

    V const b0 = V::make(0.1, 0.1,  2.0);
    V const b1 = V::make(0.1, 0.9,  2.0);
    V const b2 = V::make(0.9, 0.1,  2.0);
    V const b3 = V::make(0.2, 0.2, -2.0);

    geometry::Tetrahedron<V> const A = geometry::make_tetrahedron(a0, a1, a2 ,a3);
    geometry::Tetrahedron<V> const B = geometry::make_tetrahedron(b0, b1, b2 ,b3);

    MyCallback callback;

    bool const test = geometry::contacts_tetrahedron_tetrahedron(A, B, callback );

    BOOST_CHECK(test);

    BOOST_CHECK_EQUAL(callback.m_contacts.size(), 6u);

    BOOST_CHECK_CLOSE( callback.m_contacts[0].m_point(0), 0.244481, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[0].m_point(1), 0.244481, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[0].m_point(2), 0.742487, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[0].m_normal(0), 0.706665, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[0].m_normal(1), 0.706665, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[0].m_normal(2), 0.0353333, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[0].m_distance, -0.319536, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[1].m_point(0), 0.262902, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[1].m_point(1), 0.262902, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[1].m_point(2), 0.00564511, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[1].m_normal(0), 0.706665, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[1].m_normal(1), 0.706665, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[1].m_normal(2), 0.0353333, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[1].m_distance, -0.319536, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[2].m_point(0), 0.0305759, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[2].m_point(1), 0.48275, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[2].m_point(2), 0.255224, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[2].m_normal(0), 0.706665, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[2].m_normal(1), 0.706665, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[2].m_normal(2), 0.0353333, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[2].m_distance, -0.319536, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[3].m_point(0), 0.063152, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[3].m_point(1), 0.463152, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[3].m_point(2), -0.0043424, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[3].m_normal(0), 0.706665, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[3].m_normal(1), 0.706665, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[3].m_normal(2), 0.0353333, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[3].m_distance, -0.319536, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[4].m_point(0), 0.48275, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[4].m_point(1), 0.0305759, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[4].m_point(2), 0.255224, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[4].m_normal(0), 0.706665, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[4].m_normal(1), 0.706665, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[4].m_normal(2), 0.0353333, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[4].m_distance, -0.319536, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[5].m_point(0), 0.463152, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[5].m_point(1), 0.063152, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[5].m_point(2), -0.0043424, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[5].m_normal(0), 0.706665, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[5].m_normal(1), 0.706665, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[5].m_normal(2), 0.0353333, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[5].m_distance, -0.319536, 0.01);
  }
  // face of A is touching face of  B
  {
    V const a0 = V::make(0.0, 0.0, 0.0);
    V const a1 = V::make(1.0, 0.0, 0.0);
    V const a2 = V::make(0.0, 1.0, 0.0);
    V const a3 = V::make(0.0, 0.0, 2.0);

    V const b0 = V::make( 0.9,  0.9,  0.1);
    V const b1 = V::make( 0.9, -0.1,  0.1);
    V const b2 = V::make(-0.1,  0.9,  0.1);
    V const b3 = V::make( 0.9,  0.9, -2.0 );

    geometry::Tetrahedron<V> const A = geometry::make_tetrahedron(a0, a1, a2 ,a3);
    geometry::Tetrahedron<V> const B = geometry::make_tetrahedron(b0, b1, b2 ,b3);

    MyCallback callback;

    bool const test = geometry::contacts_tetrahedron_tetrahedron(A, B, callback );

    BOOST_CHECK(test);

    BOOST_CHECK_EQUAL(callback.m_contacts.size(), 10u);

    BOOST_CHECK_CLOSE( callback.m_contacts[0].m_point(0),   0.9,       0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[0].m_point(1),   0.05,      0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[0].m_point(2),   0.05,      0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[0].m_normal(0),  0.0,       0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[0].m_normal(1),  0.0,       0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[0].m_normal(2), -1.0,       0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[0].m_distance,  -0.1,       0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[1].m_point(0),   0.847619,  0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[1].m_point(1),   0.0,       0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[1].m_point(2),   0.05,      0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[1].m_normal(0),  0.0,       0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[1].m_normal(1),  0.0,       0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[1].m_normal(2), -1.0,       0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[1].m_distance,  -0.1,       0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[2].m_point(0),   0.9,       0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[2].m_point(1),   0.0,       0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[2].m_point(2),   0.05,      0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[2].m_normal(0),  0.0,       0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[2].m_normal(1),  0.0,       0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[2].m_normal(2), -1.0,       0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[2].m_distance,  -0.1,       0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[3].m_point(0),   0.05,      0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[3].m_point(1),   0.9,       0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[3].m_point(2),   0.05,      0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[3].m_normal(0),  0.0,       0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[3].m_normal(1),  0.0,       0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[3].m_normal(2), -1.0,       0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[3].m_distance,  -0.1,       0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[4].m_point(0),   0.0,       0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[4].m_point(1),   0.847619,  0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[4].m_point(2),   0.05,      0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[4].m_normal(0),  0.0,       0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[4].m_normal(1),  0.0,       0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[4].m_normal(2), -1.0,       0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[4].m_distance,  -0.1,       0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[5].m_point(0),   0.0,       0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[5].m_point(1),   0.9,       0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[5].m_point(2),   0.05,      0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[5].m_normal(0),  0.0,       0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[5].m_normal(1),  0.0,       0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[5].m_normal(2), -1.0,       0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[5].m_distance,  -0.1,       0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[6].m_point(0),   0.0,       0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[6].m_point(1),   0.8,       0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[6].m_point(2),   0.05,      0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[6].m_normal(0),  0.0,       0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[6].m_normal(1),  0.0,       0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[6].m_normal(2), -1.0,       0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[6].m_distance,  -0.1,       0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[7].m_point(0),   0.1,       0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[7].m_point(1),   0.9,       0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[7].m_point(2),   0.05,      0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[7].m_normal(0),  0.0,       0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[7].m_normal(1),  0.0,       0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[7].m_normal(2), -1.0,       0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[7].m_distance,  -0.1,       0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[8].m_point(0),   0.8,       0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[8].m_point(1),   0,         0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[8].m_point(2),   0.05,      0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[8].m_normal(0),  0.0,       0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[8].m_normal(1),  0.0,       0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[8].m_normal(2), -1.0,       0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[8].m_distance,  -0.1,       0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[9].m_point(0),   0.9,       0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[9].m_point(1),   0.1,       0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[9].m_point(2),   0.05,      0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[9].m_normal(0),  0.0,       0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[9].m_normal(1),  0.0,       0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[9].m_normal(2), -1.0,       0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[9].m_distance,  -0.1,       0.01);

  }

  // face of A is touching face of  B, only two surface faces
  {
    V const a0 = V::make(0.0, 0.0, 0.0);
    V const a1 = V::make(1.0, 0.0, 0.0);
    V const a2 = V::make(0.0, 1.0, 0.0);
    V const a3 = V::make(0.0, 0.0, 2.0);

    V const b0 = V::make( 0.9,  0.9,  0.1);
    V const b1 = V::make( 0.9, -0.1,  0.1);
    V const b2 = V::make(-0.1,  0.9,  0.1);
    V const b3 = V::make( 0.9,  0.9, -2.0 );

    std::vector<bool> surface_A(4u,false);
    std::vector<bool> surface_B(4u,false);

    surface_A[3] = true;
    surface_B[3] = true;

    geometry::Tetrahedron<V> const A = geometry::make_tetrahedron(a0, a1, a2 ,a3);
    geometry::Tetrahedron<V> const B = geometry::make_tetrahedron(b0, b1, b2 ,b3);

    MyCallback callback;

    bool const test = geometry::contacts_tetrahedron_tetrahedron(A, B, callback, surface_A, surface_B, geometry::MOST_OPPOSING_SURFACES() );

    BOOST_CHECK(test);

    BOOST_CHECK_EQUAL(callback.m_contacts.size(), 10u);

    BOOST_CHECK_CLOSE( callback.m_contacts[0].m_point(0),   0.9,       0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[0].m_point(1),   0.05,      0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[0].m_point(2),   0.05,      0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[0].m_normal(0),  0.0,       0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[0].m_normal(1),  0.0,       0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[0].m_normal(2), -1.0,       0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[0].m_distance,  -0.1,       0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[1].m_point(0),   0.847619,  0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[1].m_point(1),   0.0,       0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[1].m_point(2),   0.05,      0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[1].m_normal(0),  0.0,       0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[1].m_normal(1),  0.0,       0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[1].m_normal(2), -1.0,       0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[1].m_distance,  -0.1,       0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[2].m_point(0),   0.9,       0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[2].m_point(1),   0.0,       0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[2].m_point(2),   0.05,      0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[2].m_normal(0),  0.0,       0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[2].m_normal(1),  0.0,       0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[2].m_normal(2), -1.0,       0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[2].m_distance,  -0.1,       0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[3].m_point(0),   0.05,      0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[3].m_point(1),   0.9,       0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[3].m_point(2),   0.05,      0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[3].m_normal(0),  0.0,       0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[3].m_normal(1),  0.0,       0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[3].m_normal(2), -1.0,       0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[3].m_distance,  -0.1,       0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[4].m_point(0),   0.0,       0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[4].m_point(1),   0.847619,  0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[4].m_point(2),   0.05,      0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[4].m_normal(0),  0.0,       0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[4].m_normal(1),  0.0,       0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[4].m_normal(2), -1.0,       0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[4].m_distance,  -0.1,       0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[5].m_point(0),   0.0,       0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[5].m_point(1),   0.9,       0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[5].m_point(2),   0.05,      0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[5].m_normal(0),  0.0,       0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[5].m_normal(1),  0.0,       0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[5].m_normal(2), -1.0,       0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[5].m_distance,  -0.1,       0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[6].m_point(0),   0.0,       0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[6].m_point(1),   0.8,       0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[6].m_point(2),   0.05,      0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[6].m_normal(0),  0.0,       0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[6].m_normal(1),  0.0,       0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[6].m_normal(2), -1.0,       0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[6].m_distance,  -0.1,       0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[7].m_point(0),   0.1,       0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[7].m_point(1),   0.9,       0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[7].m_point(2),   0.05,      0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[7].m_normal(0),  0.0,       0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[7].m_normal(1),  0.0,       0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[7].m_normal(2), -1.0,       0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[7].m_distance,  -0.1,       0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[8].m_point(0),   0.8,       0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[8].m_point(1),   0,         0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[8].m_point(2),   0.05,      0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[8].m_normal(0),  0.0,       0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[8].m_normal(1),  0.0,       0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[8].m_normal(2), -1.0,       0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[8].m_distance,  -0.1,       0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[9].m_point(0),   0.9,       0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[9].m_point(1),   0.1,       0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[9].m_point(2),   0.05,      0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[9].m_normal(0),  0.0,       0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[9].m_normal(1),  0.0,       0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[9].m_normal(2), -1.0,       0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[9].m_distance,  -0.1,       0.01);
    
  }

  // Edge from B goes through A and vice-versa
  {
    V const a0 = V::make(0.0, 0.0, 0.0);
    V const a1 = V::make(1.0, 0.0, 0.0);
    V const a2 = V::make(0.0, 1.0, 0.0);
    V const a3 = V::make(0.0, 0.0, 1.0);

    V const b0 = V::make(0.3, 0.3, -1.0);
    V const b1 = V::make(1.3, 0.3,  0.0);
    V const b2 = V::make(0.3, 1.3,  0.0);
    V const b3 = V::make(0.3, 0.3,  1.0);

    geometry::Tetrahedron<V> const A = geometry::make_tetrahedron(a0, a1, a2 ,a3);
    geometry::Tetrahedron<V> const B = geometry::make_tetrahedron(b0, b1, b2 ,b3);

    MyCallback callback;

    bool const test = geometry::contacts_tetrahedron_tetrahedron(A, B, callback );

    BOOST_CHECK(test);

    BOOST_CHECK_EQUAL(callback.m_contacts.size(), 4u);

    BOOST_CHECK_CLOSE( callback.m_contacts[0].m_point(0), 0.4, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[0].m_point(1), 0.4, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[0].m_point(2), 0.4, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[0].m_normal(0), 0.707107, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[0].m_normal(1), 0.707107, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[0].m_normal(2), -0, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[0].m_distance, -0.282843, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[1].m_point(0), 0.4, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[1].m_point(1), 0.4, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[1].m_point(2), 0, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[1].m_normal(0), 0.707107, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[1].m_normal(1), 0.707107, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[1].m_normal(2), -0, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[1].m_distance, -0.282843, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[2].m_point(0), 0.2, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[2].m_point(1), 0.6, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[2].m_point(2), 0, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[2].m_normal(0), 0.707107, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[2].m_normal(1), 0.707107, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[2].m_normal(2), -0, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[2].m_distance, -0.282843, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[3].m_point(0), 0.6, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[3].m_point(1), 0.2, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[3].m_point(2), 0, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[3].m_normal(0), 0.707107, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[3].m_normal(1), 0.707107, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[3].m_normal(2), -0, 0.01);
    BOOST_CHECK_CLOSE( callback.m_contacts[3].m_distance, -0.282843, 0.01);
    
  }

}

BOOST_AUTO_TEST_SUITE_END();
