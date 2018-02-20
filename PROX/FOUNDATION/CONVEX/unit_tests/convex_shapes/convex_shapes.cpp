#include <tiny.h>
#include <geometry.h>   // needed for geometry::Box, geometry::Sphere
#include <convex_shapes.h>

#define BOOST_AUTO_TEST_MAIN
#include <boost/test/auto_unit_test.hpp>
#include <boost/test/unit_test_suite.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/test_tools.hpp>

BOOST_AUTO_TEST_SUITE(convex_shapes);

BOOST_AUTO_TEST_CASE(case_by_case_testing)
{
  typedef tiny::MathTypes<double> M;
  typedef M::vector3_type         V;

  geometry::Box<V>        const box;
  convex::Capsule<M>      const capsule;
  convex::Cone<M>         const cone;
  convex::ConvexHull<M>         hull;
  convex::Cylinder<M>     const cylinder;
  convex::Ellipsoid<M>    const ellipsoid;
  geometry::Sphere<V>     const sphere;

  V const zero     = V::make( 0.0, 0.0, 0.0);
  V const up       = V::make( 0.0, 0.0, 1.0);
  V const down     = V::make( 0.0, 0.0,-1.0);
  V const right    = V::make( 1.0, 0.0, 0.0);
  V const left     = V::make(-1.0, 0.0, 0.0);
  V const forward  = V::make( 0.0, 1.0, 0.0);
  V const backward = V::make( 0.0,-1.0, 0.0);

  // Test what happens if zero vector is given as search direction
  {
    V p = box.get_support_point( zero );
    BOOST_CHECK( (p(0)!=0.0) || (p(1)!=0) || (p(2)!=0) );

    p = capsule.get_support_point( zero );
    BOOST_CHECK( (p(0)!=0.0) || (p(1)!=0) || (p(2)!=0) );

    p = cone.get_support_point( zero );
    BOOST_CHECK( (p(0)!=0.0) || (p(1)!=0) || (p(2)!=0) );

    p = cylinder.get_support_point( zero );
    BOOST_CHECK( (p(0)!=0.0) || (p(1)!=0) || (p(2)!=0) );

    p = ellipsoid.get_support_point( zero );
    BOOST_CHECK( (p(0)!=0.0) || (p(1)!=0) || (p(2)!=0) );

    p = sphere.get_support_point( zero );
    BOOST_CHECK( (p(0)!=0.0) || (p(1)!=0) || (p(2)!=0) );
  }

  // Test what happens if up vector is given
  {
    V p = box.get_support_point( up );

    BOOST_CHECK( p(0) <=  box.half_extent()(0) );
    BOOST_CHECK( p(1) <=  box.half_extent()(1) );
    BOOST_CHECK( p(0) >= -box.half_extent()(0) );
    BOOST_CHECK( p(1) >= -box.half_extent()(1) );
    BOOST_CHECK( p(2) ==  box.half_extent()(2) );
  }
  {
    V p = capsule.get_support_point( up );
    BOOST_CHECK( p(0) == 0.0 );
    BOOST_CHECK( p(1) == 0.0 );
    BOOST_CHECK( p(2) == 2.0 );
  }
  {
    V p = cone.get_support_point( up );
    BOOST_CHECK( p(0) == 0.0 );
    BOOST_CHECK( p(1) == 0.0 );
    BOOST_CHECK( p(2) == 1.0 );
  }
  {
    V p = cylinder.get_support_point( up );
    BOOST_CHECK( p(0) == 0.0 );
    BOOST_CHECK( p(1) == 0.0 );
    BOOST_CHECK( p(2) == 1.0 );
  }
  {
    V p = ellipsoid.get_support_point( up );
    BOOST_CHECK( p(0) == 0.0 );
    BOOST_CHECK( p(1) == 0.0 );
    BOOST_CHECK( p(2) == 1.0 );
  }
  {
    V p = sphere.get_support_point( up );
    BOOST_CHECK( p(0) == 0.0 );
    BOOST_CHECK( p(1) == 0.0 );
    BOOST_CHECK( p(2) == 1.0 );
  }

  // Test what happen if right vector is given
  {
    V p = box.get_support_point( right );

    BOOST_CHECK( p(0) ==  box.half_extent()(0) );
    BOOST_CHECK( p(1) <=  box.half_extent()(1) );
    BOOST_CHECK( p(2) <=  box.half_extent()(2) );
    BOOST_CHECK( p(1) >= -box.half_extent()(1) );
    BOOST_CHECK( p(2) >= -box.half_extent()(2) );
  }
  {
    V p = capsule.get_support_point( right );
    BOOST_CHECK( p(0) == 1.0 );
    BOOST_CHECK( p(1) == 0.0 );
    BOOST_CHECK( p(2) == 0.0 );
  }
  {
    V p = cone.get_support_point( right );
    BOOST_CHECK( p(0) ==  1.0 );
    BOOST_CHECK( p(1) ==  0.0 );
    BOOST_CHECK( p(2) == -1.0 );
  }
  {
    V p = cylinder.get_support_point( right );
    BOOST_CHECK( p(0) == cylinder.radius() );
    BOOST_CHECK( p(1) == 0.0 );
    BOOST_CHECK( p(2) <=  cylinder.half_height() );
    BOOST_CHECK( p(2) >= -cylinder.half_height() );
  }
  {
    V p = ellipsoid.get_support_point( right );
    BOOST_CHECK( p(0) ==  1.0 );
    BOOST_CHECK( p(1) ==  0.0 );
    BOOST_CHECK( p(2) ==  0.0 );
  }
  {
    V p = sphere.get_support_point( right );
    BOOST_CHECK( p(0) ==  1.0 );
    BOOST_CHECK( p(1) ==  0.0 );
    BOOST_CHECK( p(2) ==  0.0 );
  }

  // Test what happen if down vector is given
  {
    V p = box.get_support_point( down );

    BOOST_CHECK( p(2) == -box.half_extent()(2) );
    BOOST_CHECK( p(0) <=  box.half_extent()(0) );
    BOOST_CHECK( p(1) <=  box.half_extent()(1) );
    BOOST_CHECK( p(0) >= -box.half_extent()(0) );
    BOOST_CHECK( p(1) >= -box.half_extent()(1) );
  }
  {
    V p = capsule.get_support_point( down );
    BOOST_CHECK( p(0) ==   0.0 );
    BOOST_CHECK( p(1) ==   0.0 );
    BOOST_CHECK( p(2) ==  -2.0 );
  }
  {
    V p = cone.get_support_point( down );
    BOOST_CHECK( p(0) ==   0.0 );
    BOOST_CHECK( p(1) ==   0.0 );
    BOOST_CHECK( p(2) ==  -1.0 );
  }
  {
    V p = cylinder.get_support_point( down );
    BOOST_CHECK( p(0) == 0.0 );
    BOOST_CHECK( p(1) == 0.0 );
    BOOST_CHECK( p(2) == -cylinder.half_height() );
  }
  {
    V p = ellipsoid.get_support_point( down );
    BOOST_CHECK( p(0) ==   0.0 );
    BOOST_CHECK( p(1) ==   0.0 );
    BOOST_CHECK( p(2) ==  -1.0 );
  }
  {
    V p = sphere.get_support_point( down );
    BOOST_CHECK( p(0) ==   0.0 );
    BOOST_CHECK( p(1) ==   0.0 );
    BOOST_CHECK( p(2) ==  -1.0 );
  }
  // Test the convex hull
  {
    hull.add_point( left     );
    hull.add_point( right    );
    hull.add_point( backward );
    hull.add_point( forward  );
    hull.add_point( down     );
    hull.add_point( up       );

    V p = hull.get_support_point(right );
    BOOST_CHECK( p(0) ==   1.0 );
    BOOST_CHECK( p(1) ==   0.0 );
    BOOST_CHECK( p(2) ==   0.0 );

    p = hull.get_support_point( left );
    BOOST_CHECK( p(0) ==  -1.0 );
    BOOST_CHECK( p(1) ==   0.0 );
    BOOST_CHECK( p(2) ==   0.0 );

    p = hull.get_support_point( forward );
    BOOST_CHECK( p(0) ==   0.0 );
    BOOST_CHECK( p(1) ==   1.0 );
    BOOST_CHECK( p(2) ==   0.0 );

    p = hull.get_support_point( backward );
    BOOST_CHECK( p(0) ==   0.0 );
    BOOST_CHECK( p(1) ==  -1.0 );
    BOOST_CHECK( p(2) ==   0.0 );

    p = hull.get_support_point( up );
    BOOST_CHECK( p(0) ==   0.0 );
    BOOST_CHECK( p(1) ==   0.0 );
    BOOST_CHECK( p(2) ==   1.0 );

    p = hull.get_support_point( down );
    BOOST_CHECK( p(0) ==   0.0 );
    BOOST_CHECK( p(1) ==   0.0 );
    BOOST_CHECK( p(2) ==  -1.0 );
  }
}

BOOST_AUTO_TEST_SUITE_END();
