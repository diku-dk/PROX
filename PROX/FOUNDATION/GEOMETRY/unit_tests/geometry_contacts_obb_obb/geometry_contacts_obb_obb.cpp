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
typedef MT::value_traits       VT;

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

BOOST_AUTO_TEST_CASE(contacts_obb_obb_test)
{
  // touching bottom-top faces
  {
    V const centerA   = V::make(0.0, 0.0, 0.0);
    V const half_extA = V::make(1.0, 1.0, 1.0);
    Q const qA        = Q::identity();

    geometry::OBB<MT> obbA = geometry::make_obb<MT>(centerA, qA, half_extA);

    V const centerB   = V::make(0.0, -2.0, 0.0);
    V const half_extB = V::make(1.0, 1.0, 1.0);
    Q const qB        = Q::Ry(VT::pi_quarter());

    geometry::OBB<MT> obbB = geometry::make_obb<MT>(centerB, qB, half_extB);

    MyCallback callback1;
    MyCallback callback2;

    bool const test1 = geometry::contacts_obb_obb(obbA,obbB, 0.0, callback1);
    bool const test2 = geometry::contacts_obb_obb(obbB,obbA, 0.0, callback2);

    BOOST_CHECK(test1);
    BOOST_CHECK(test2);

    BOOST_CHECK_EQUAL(callback1.m_contacts.size(), 8u);
    BOOST_CHECK_EQUAL(callback2.m_contacts.size(), 8u);

    for(unsigned int k = 0u;k <8u;++k)
    {
      BOOST_CHECK( fabs(       callback1.m_contacts[k].m_normal(2)) < 0.00001 );
      BOOST_CHECK( fabs(       callback1.m_contacts[k].m_normal(0)) < 0.00001 );
      BOOST_CHECK( 1.0  - fabs(callback1.m_contacts[k].m_normal(1)) < 0.00001 );
      BOOST_CHECK(       fabs( callback1.m_contacts[k].m_distance)  < 0.00001 );
    }


  }
  // touching left-right faces
  {
    V const centerA   = V::make(0.0, 0.0, 0.0);
    V const half_extA = V::make(1.0, 1.0, 1.0);
    Q const qA        = Q::identity();

    geometry::OBB<MT> obbA = geometry::make_obb<MT>(centerA, qA, half_extA);

    V const centerB   = V::make(2.0, 0.0, 0.0);
    V const half_extB = V::make(1.0, 1.0, 1.0);
    Q const qB        = Q::Rx(VT::pi_quarter());

    geometry::OBB<MT> obbB = geometry::make_obb<MT>(centerB, qB, half_extB);

    MyCallback callback1;
    MyCallback callback2;

    bool const test1 = geometry::contacts_obb_obb(obbA,obbB, 0.0, callback1);
    bool const test2 = geometry::contacts_obb_obb(obbB,obbA, 0.0, callback2);

    BOOST_CHECK(test1);
    BOOST_CHECK(test2);

    BOOST_CHECK_EQUAL(callback1.m_contacts.size(), 8u);
    BOOST_CHECK_EQUAL(callback2.m_contacts.size(), 8u);

    for(unsigned int k = 0u;k <8u;++k)
    {
      BOOST_CHECK( fabs(       callback1.m_contacts[k].m_normal(2)) < 0.00001 );
      BOOST_CHECK( fabs(       callback1.m_contacts[k].m_normal(1)) < 0.00001 );
      BOOST_CHECK( 1.0  - fabs(callback1.m_contacts[k].m_normal(0)) < 0.00001 );
      BOOST_CHECK(       fabs( callback1.m_contacts[k].m_distance)  < 0.00001 );
    }


  }
  // touching front-back faces
  {
    V const centerA   = V::make(0.0, 0.0, 0.0);
    V const half_extA = V::make(1.0, 1.0, 1.0);
    Q const qA        = Q::identity();

    geometry::OBB<MT> obbA = geometry::make_obb<MT>(centerA, qA, half_extA);

    V const centerB   = V::make(0.0, 0.0, 2.0);
    V const half_extB = V::make(1.0, 1.0, 1.0);
    Q const qB        = Q::Rz(VT::pi_quarter());

    geometry::OBB<MT> obbB = geometry::make_obb<MT>(centerB, qB, half_extB);

    MyCallback callback1;
    MyCallback callback2;

    bool const test1 = geometry::contacts_obb_obb(obbA,obbB, 0.0, callback1);
    bool const test2 = geometry::contacts_obb_obb(obbB,obbA, 0.0, callback2);

    BOOST_CHECK(test1);
    BOOST_CHECK(test2);

    BOOST_CHECK_EQUAL(callback1.m_contacts.size(), 8u);
    BOOST_CHECK_EQUAL(callback2.m_contacts.size(), 8u);

    for(unsigned int k = 0u;k <8u;++k)
    {
      BOOST_CHECK( fabs(       callback1.m_contacts[k].m_normal(0)) < 0.00001 );
      BOOST_CHECK( fabs(       callback1.m_contacts[k].m_normal(1)) < 0.00001 );
      BOOST_CHECK( 1.0  - fabs(callback1.m_contacts[k].m_normal(2)) < 0.00001 );
      BOOST_CHECK(       fabs( callback1.m_contacts[k].m_distance)  < 0.00001 );
    }


  }
  // separating bottom-top faces
  {
    V const centerA   = V::make(0.0, 0.0, 0.0);
    V const half_extA = V::make(1.0, 1.0, 1.0);
    Q const qA        = Q::identity();

    geometry::OBB<MT> obbA = geometry::make_obb<MT>(centerA, qA, half_extA);

    V const centerB   = V::make(0.0, -2.1, 0.0);
    V const half_extB = V::make(1.0, 1.0, 1.0);
    Q const qB        = Q::Ry(VT::pi_quarter());

    geometry::OBB<MT> obbB = geometry::make_obb<MT>(centerB, qB, half_extB);

    MyCallback callback1;
    MyCallback callback2;

    bool const test1 = geometry::contacts_obb_obb(obbA,obbB, 0.0, callback1);
    bool const test2 = geometry::contacts_obb_obb(obbB,obbA, 0.0, callback2);

    BOOST_CHECK(!test1);
    BOOST_CHECK(!test2);

    BOOST_CHECK_EQUAL(callback1.m_contacts.size(), 0u);
    BOOST_CHECK_EQUAL(callback2.m_contacts.size(), 0u);

  }
  // separating left-right faces
  {
    V const centerA   = V::make(0.0, 0.0, 0.0);
    V const half_extA = V::make(1.0, 1.0, 1.0);
    Q const qA        = Q::identity();

    geometry::OBB<MT> obbA = geometry::make_obb<MT>(centerA, qA, half_extA);

    V const centerB   = V::make(2.1, 0.0, 0.0);
    V const half_extB = V::make(1.0, 1.0, 1.0);
    Q const qB        = Q::Rx(VT::pi_quarter());

    geometry::OBB<MT> obbB = geometry::make_obb<MT>(centerB, qB, half_extB);

    MyCallback callback1;
    MyCallback callback2;

    bool const test1 = geometry::contacts_obb_obb(obbA,obbB, 0.0, callback1);
    bool const test2 = geometry::contacts_obb_obb(obbB,obbA, 0.0, callback2);

    BOOST_CHECK(!test1);
    BOOST_CHECK(!test2);

    BOOST_CHECK_EQUAL(callback1.m_contacts.size(), 0u);
    BOOST_CHECK_EQUAL(callback2.m_contacts.size(), 0u);

  }
  // separating front-back faces
  {
    V const centerA   = V::make(0.0, 0.0, 0.0);
    V const half_extA = V::make(1.0, 1.0, 1.0);
    Q const qA        = Q::identity();

    geometry::OBB<MT> obbA = geometry::make_obb<MT>(centerA, qA, half_extA);

    V const centerB   = V::make(0.0, 0.0, 2.1);
    V const half_extB = V::make(1.0, 1.0, 1.0);
    Q const qB        = Q::Rz(VT::pi_quarter());

    geometry::OBB<MT> obbB = geometry::make_obb<MT>(centerB, qB, half_extB);

    MyCallback callback1;
    MyCallback callback2;

    bool const test1 = geometry::contacts_obb_obb(obbA,obbB, 0.0, callback1);
    bool const test2 = geometry::contacts_obb_obb(obbB,obbA, 0.0, callback2);

    BOOST_CHECK(!test1);
    BOOST_CHECK(!test2);

    BOOST_CHECK_EQUAL(callback1.m_contacts.size(), 0u);
    BOOST_CHECK_EQUAL(callback2.m_contacts.size(), 0u);

  }
  // A inside of B
  {
    V const centerA   = V::make(0.0,   0.5,  0.0 );
    V const half_extA = V::make(0.5,   0.25, 0.5);
    Q const qA        = Q::identity();

    geometry::OBB<MT> obbA = geometry::make_obb<MT>(centerA, qA, half_extA);

    V const centerB   = V::make(0.0, 0.0, 0.0);
    V const half_extB = V::make(1.0, 1.0, 1.0);
    Q const qB        = Q::identity();

    geometry::OBB<MT> obbB = geometry::make_obb<MT>(centerB, qB, half_extB);

    MyCallback callback1;
    MyCallback callback2;

    bool const test1 = geometry::contacts_obb_obb(obbA,obbB, 0.0, callback1);
    bool const test2 = geometry::contacts_obb_obb(obbB,obbA, 0.0, callback2);

    BOOST_CHECK(test1);
    BOOST_CHECK(test2);

    BOOST_CHECK_EQUAL(callback1.m_contacts.size(), 4u);
    BOOST_CHECK_EQUAL(callback2.m_contacts.size(), 4u);

    for(unsigned int k = 0u;k <4u;++k)
    {
      BOOST_CHECK( 0.5 - fabs(callback1.m_contacts[k].m_point(0))   < 0.00001 );
      BOOST_CHECK( 0.5 - fabs(callback1.m_contacts[k].m_point(1))   < 0.00001 );
      BOOST_CHECK( 0.5 - fabs(callback1.m_contacts[k].m_point(2))   < 0.00001 );
      BOOST_CHECK( fabs(       callback1.m_contacts[k].m_normal(0)) < 0.00001 );
      BOOST_CHECK( 1.0  - fabs(callback1.m_contacts[k].m_normal(1)) < 0.00001 );
      BOOST_CHECK( fabs(       callback1.m_contacts[k].m_normal(2)) < 0.00001 );
      BOOST_CHECK( 0.5 - fabs( callback1.m_contacts[k].m_distance)  < 0.00001 ); // depth is off due to technicality in method, should have been 0.75
    }

  }
  // separating edge-edge-case
  {
    V const centerA   = V::make(0.0, 0.0, 0.0);
    V const half_extA = V::make(1.0, 1.0, 1.0);
    Q const qA        = Q::identity();

    geometry::OBB<MT> obbA = geometry::make_obb<MT>(centerA, qA, half_extA);

    V const centerB   = V::make(2.01, 2.01, 0.0);
    V const half_extB = V::make(1.0, 1.0, 1.0);
    Q const qB        = Q::Ru(VT::pi_half(), centerB );

    geometry::OBB<MT> obbB = geometry::make_obb<MT>(centerB, qB, half_extB);

    MyCallback callback1;
    MyCallback callback2;

    bool const test1 = geometry::contacts_obb_obb(obbA,obbB, 0.0, callback1);
    bool const test2 = geometry::contacts_obb_obb(obbB,obbA, 0.0, callback2);

    BOOST_CHECK(!test1);
    BOOST_CHECK(!test2);

    BOOST_CHECK_EQUAL(callback1.m_contacts.size(), 0u);
    BOOST_CHECK_EQUAL(callback2.m_contacts.size(), 0u);
  }
  // touching edge-edge-case
  {
    V const centerA   = V::make(0.0, 0.0, 0.0);
    V const half_extA = V::make(1.0, 1.0, 1.0);
    Q const qA        = Q::identity();

    geometry::OBB<MT> obbA = geometry::make_obb<MT>(centerA, qA, half_extA);


    V const centerB   = V::make(1.999999, 1.999999, 0.0);
    V const half_extB = V::make(1.0, 1.0, 1.0);
    Q const qB        = Q::Ru(VT::pi_half(), centerB );

    geometry::OBB<MT> obbB = geometry::make_obb<MT>(centerB, qB, half_extB);

    MyCallback callback1;
    MyCallback callback2;

    bool const test1 = geometry::contacts_obb_obb(obbA,obbB, 0.0, callback1);
    bool const test2 = geometry::contacts_obb_obb(obbB,obbA, 0.0, callback2);

    BOOST_CHECK(test1);
    BOOST_CHECK(test2);

    BOOST_CHECK_EQUAL(callback1.m_contacts.size(), 4u);
    BOOST_CHECK_EQUAL(callback2.m_contacts.size(), 4u);


    BOOST_CHECK( fabs(1.0 - callback1.m_contacts[0].m_point(0)) < 0.00001);
    BOOST_CHECK( fabs(1.0 - callback1.m_contacts[0].m_point(1)) < 0.00001);
    BOOST_CHECK( fabs(      callback1.m_contacts[0].m_point(2)) < 0.00001);
    BOOST_CHECK( fabs(sqrt(0.5) - callback1.m_contacts[0].m_normal(0)) < 0.00001);
    BOOST_CHECK( fabs(sqrt(0.5) - callback1.m_contacts[0].m_normal(1)) < 0.00001);
    BOOST_CHECK( fabs(            callback1.m_contacts[0].m_normal(2)) < 0.00001);
    BOOST_CHECK( fabs( callback1.m_contacts[0].m_distance) < 0.00001);


    BOOST_CHECK( fabs(1.0 - callback1.m_contacts[1].m_point(0)) < 0.00001);
    BOOST_CHECK( fabs(1.0 - callback1.m_contacts[1].m_point(1)) < 0.00001);
    BOOST_CHECK( fabs(      callback1.m_contacts[1].m_point(2)) < 0.00001);
    BOOST_CHECK( fabs(sqrt(0.5) - callback1.m_contacts[1].m_normal(0)) < 0.00001);
    BOOST_CHECK( fabs(sqrt(0.5) - callback1.m_contacts[1].m_normal(1)) < 0.00001);
    BOOST_CHECK( fabs(            callback1.m_contacts[1].m_normal(2)) < 0.00001);
    BOOST_CHECK( fabs( callback1.m_contacts[1].m_distance) < 0.00001);

    BOOST_CHECK( fabs(1.0 - callback1.m_contacts[2].m_point(0)) < 0.00001);
    BOOST_CHECK( fabs(1.0 - callback1.m_contacts[2].m_point(1)) < 0.00001);
    BOOST_CHECK( fabs(      callback1.m_contacts[2].m_point(2)) < 0.00001);
    BOOST_CHECK( fabs(sqrt(0.5) - callback1.m_contacts[2].m_normal(0)) < 0.00001);
    BOOST_CHECK( fabs(sqrt(0.5) - callback1.m_contacts[2].m_normal(1)) < 0.00001);
    BOOST_CHECK( fabs(            callback1.m_contacts[2].m_normal(2)) < 0.00001);
    BOOST_CHECK( fabs( callback1.m_contacts[2].m_distance) < 0.00001);

    BOOST_CHECK( fabs(1.0 - callback1.m_contacts[3].m_point(0)) < 0.00001);
    BOOST_CHECK( fabs(1.0 - callback1.m_contacts[3].m_point(1)) < 0.00001);
    BOOST_CHECK( fabs(      callback1.m_contacts[3].m_point(2)) < 0.00001);
    BOOST_CHECK( fabs(sqrt(0.5) - callback1.m_contacts[3].m_normal(0)) < 0.00001);
    BOOST_CHECK( fabs(sqrt(0.5) - callback1.m_contacts[3].m_normal(1)) < 0.00001);
    BOOST_CHECK( fabs(            callback1.m_contacts[3].m_normal(2)) < 0.00001);
    BOOST_CHECK( fabs( callback1.m_contacts[3].m_distance) < 0.00001);

  }

}

BOOST_AUTO_TEST_SUITE_END();
