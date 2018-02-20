#include <tiny.h>
#include <geometry.h>
#include <convex.h>

#define BOOST_AUTO_TEST_MAIN
#include <boost/test/auto_unit_test.hpp>
#include <boost/test/unit_test_suite.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/test_tools.hpp>

#include <cmath>

typedef tiny::MathTypes<double>  M;
typedef M::quaternion_type       Q;
typedef M::vector3_type          V;
typedef M::real_type             T;
typedef M::coordsys_type         C;
typedef M::value_traits          VT;

BOOST_AUTO_TEST_SUITE(convex_get_scale);

BOOST_AUTO_TEST_CASE(get_scale)
{
  geometry::Box<V>          box;
  convex::Capsule<M>          capsule;
  convex::Cone<M>             cone;
  convex::ConvexHull<M>       hull;
  convex::Cylinder<M>         cylinder;
  convex::Ellipsoid<M>        ellipsoid;
  geometry::Sphere<V>         sphere;

  //get_scale is supposed to return the minimum extent of the shape
  //spheres
  {
    sphere.radius() = 0.75f;
    BOOST_CHECK( sphere.get_scale() == 1.5f); //for a sphere, the scale is the diameter
    sphere.radius() = 0.5f;
    BOOST_CHECK( sphere.get_scale() == 1.0f); //for a sphere, the scale is the diameter
  }
  //boxes
  {
    box.half_extent()(0) = 0.5f;//
    box.half_extent()(1) = 0.5f;// a unit cube, scale should be 1.0f
    box.half_extent()(2) = 0.5f;//
    BOOST_CHECK( box.get_scale() == 1.0f);
    box.half_extent()(1) = 0.2f; //decreasing the extent should change the scale
    BOOST_CHECK(box.get_scale() == 0.40f);
    box.half_extent()(2) = 1.2f; //this should not change the scale
    BOOST_CHECK( box.get_scale() == 0.40f);
    box.half_extent()(2) = 0.0f; //dropping one dimension, this should not change the scale
    BOOST_CHECK( box.get_scale() == 0.40f);
  }
  //capsule
  {
    capsule.half_height() = 1.50f;
    capsule.radius() = 0.7f;
    BOOST_CHECK( capsule.get_scale() == 1.4f);
    capsule.radius() = 0.5f;
    BOOST_CHECK( capsule.get_scale() == 1.0f);
    capsule.half_height() = 0.25f;
    BOOST_CHECK( capsule.get_scale() == 0.5f);
    capsule.radius() = 2.5f;
    BOOST_CHECK( capsule.get_scale() == 0.5f);
  }
  //cone
  {
    cone.half_height() = 1.50f;
    cone.base_radius() = 0.7f;
    BOOST_CHECK( cone.get_scale() == 1.4f);
    cone.base_radius() = 0.5f;
    BOOST_CHECK( cone.get_scale() == 1.0f);
    cone.half_height() = 0.25f;
    BOOST_CHECK( cone.get_scale() == 0.5f);
    cone.base_radius() = 2.5f;
    BOOST_CHECK( cone.get_scale() == 0.5f);
    cone.base_radius() = 2.5f;
    BOOST_CHECK( cone.get_scale() == 0.5f);
  }
  //cylinder
  {
    cylinder.half_height() = 1.50f;
    cylinder.radius() = 0.7f;
    BOOST_CHECK( cylinder.get_scale() == 1.4f);
    cylinder.radius() = 0.5f;
    BOOST_CHECK( cylinder.get_scale() == 1.0f);
    cylinder.half_height() = 0.25f;
    BOOST_CHECK( cylinder.get_scale() == 0.5f);
    cylinder.radius() = 2.5f;
    BOOST_CHECK( cylinder.get_scale() == 0.5f);
  }
  //ellipsoid
  {
    ellipsoid.scale()(0) = 1.40f;
    ellipsoid.scale()(1) = 0.7f;
    ellipsoid.scale()(2) = 0.35f;
    BOOST_CHECK( ellipsoid.get_scale() == 0.7f);
    ellipsoid.scale()(2) = 0.7f;
    BOOST_CHECK( ellipsoid.get_scale() == 1.4f);
    ellipsoid.scale()(2) = 0.9f;
    BOOST_CHECK( ellipsoid.get_scale() == 1.4f);
    ellipsoid.scale()(2) = 0.0f; //dropping one dimension, should not change the scale
    BOOST_CHECK( ellipsoid.get_scale() == 1.4f);
  }
  //convex hull
  {
    hull.add_point( V::make( 1.0,  4.0,  6.0 ) );
    hull.add_point( V::make( 1.0,  6.0,  8.0 ) );
    hull.add_point( V::make( 1.0,  7.0,  4.0 ) );
    hull.add_point( V::make( 1.0,  8.0,  7.0 ) );
    hull.add_point( V::make( 3.0,  5.0,  9.0 ) );
    hull.add_point( V::make( 3.0,  7.0,  5.0 ) );
    hull.add_point( V::make( 3.0,  9.0,  7.0 ) );
    hull.add_point( V::make( 4.0,  7.0,  9.0 ) );
    hull.add_point( V::make( 4.0,  9.0,  6.0 ) );
    hull.add_point( V::make( 5.0,  6.0,  9.0 ) );
    hull.add_point( V::make( 5.0,  7.0,  8.0 ) );
    hull.add_point( V::make( 5.0,  8.0,  6.0 ) );
    BOOST_CHECK(hull.get_scale() == 4.0f);
    
    hull.add_point( V::make( 10.0, 10.0, 10.0 ) );
    BOOST_CHECK(hull.get_scale() == 6.0f);
    
    hull.add_point( V::make( 1.0,  4.0,  4.0 ) );
    BOOST_CHECK(hull.get_scale() == 6.0f);
  }
}

BOOST_AUTO_TEST_SUITE_END();
