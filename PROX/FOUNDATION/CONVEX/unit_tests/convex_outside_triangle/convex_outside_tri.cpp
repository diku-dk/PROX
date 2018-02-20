#include <tiny.h>
#include <convex.h>

#define BOOST_AUTO_TEST_MAIN
#include <boost/test/auto_unit_test.hpp>
#include <boost/test/unit_test_suite.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/test_tools.hpp>

BOOST_AUTO_TEST_SUITE(convex_outside_triangle);

BOOST_AUTO_TEST_CASE(case_by_case_test)
{
  typedef tiny::MathTypes<double>       math_types;
  typedef math_types::vector3_type     V;

  V const a = V::make(0.0, 0.0, 0.0);
  V const b = V::make(1.0, 0.0, 0.0);;
  V const c = V::make(0.0, 1.0, 0.0);;
  V const q = V::make(0.33, 0.33, -1.0);

  {
    V p = V::make(0.33, 0.33,  1.0);
    bool outside = convex::outside_triangle(p, a, b, c, q );
    BOOST_CHECK( outside );
  }
  {
    V p = V::make(0.1, 0.1,  -1.0);
    bool outside = convex::outside_triangle(p, a, b, c, q );
    BOOST_CHECK( !outside );
  }
  {
    V p = V::make(0.1, 0.1,  0.0);
    bool outside = convex::outside_triangle(p, a, b, c, q );
    BOOST_CHECK( outside );
  }
}

BOOST_AUTO_TEST_SUITE_END();
