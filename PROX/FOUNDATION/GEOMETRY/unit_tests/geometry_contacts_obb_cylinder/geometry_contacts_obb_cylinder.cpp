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

BOOST_AUTO_TEST_CASE(contacts_obb_cylinder_test)
{
  {
    V const center    = V::make(0.0, 0.0, 0.0);
    V const half_ext  = V::make(1.0, 1.0, 1.0);
    Q const q         = Q::identity();

    V const center2   = V::zero();
    V const axis      = V::i();
    T const height    = 1.0;
    T const radius    = 1.0;

    geometry::OBB<MT>     const & obb = geometry::make_obb<MT>(center, q, half_ext);
    geometry::Cylinder<V> const & cyl = geometry::make_cylinder(radius, height, axis, center2);

    MyCallback mycallback;

    bool const test = geometry::contacts_obb_cylinder(obb, cyl, 0.0, mycallback, false );

    BOOST_CHECK(!test);
  }
}

BOOST_AUTO_TEST_SUITE_END();
