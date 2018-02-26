#include <hyper.h>

#define BOOST_AUTO_TEST_MAIN
#include <boost/test/auto_unit_test.hpp>
#include <boost/test/unit_test_suite.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/test_tools.hpp>


typedef hyper::MathPolicy<float> MT;
typedef MT::real_type          T;
typedef MT::vector3_type       V;
typedef MT::matrix3x3_type     M;
typedef MT::quaternion_type    Q;
typedef MT::coordsys_type      C;
typedef MT::value_traits       VT;



BOOST_AUTO_TEST_SUITE(hyper);

BOOST_AUTO_TEST_CASE(implicit_surface)
{
  // Add your test code here


  hyper::BoxSurface<MT> box = make_box<MT>(V::make(0.0, 0.0, 0.0),V::make(2.0, 2.0, 2.0));


  // center point
  {
    T const d = box( V::make( 0.0, 0.0, 0.0) );
    BOOST_CHECK_CLOSE( d, -1.0, 0.01 );
  }

  // corner points
  {
    T const d = box( V::make(  1.0,  1.0,  1.0) );
    BOOST_CHECK_CLOSE( d, 0.0, 0.01 );
  }
  {
    T const d = box( V::make(  1.0,  1.0, -1.0) );
    BOOST_CHECK_CLOSE( d, 0.0, 0.01 );
  }
  {
    T const d = box( V::make(  1.0, -1.0,  1.0) );
    BOOST_CHECK_CLOSE( d, 0.0, 0.01 );
  }
  {
    T const d = box( V::make(  1.0, -1.0, -1.0) );
    BOOST_CHECK_CLOSE( d, 0.0, 0.01 );
  }
  {
    T const d = box( V::make( -1.0,  1.0,  1.0) );
    BOOST_CHECK_CLOSE( d, 0.0, 0.01 );
  }


  // edge points
  {
    T const d = box( V::make(  1.0,  1.0,  0.0) );
    BOOST_CHECK_CLOSE( d, 0.0, 0.01 );
  }
  {
    T const d = box( V::make(  0.0,  1.0,  1.0) );
    BOOST_CHECK_CLOSE( d, 0.0, 0.01 );
  }
  {
    T const d = box( V::make( -1.0, -1.0,  0.0) );
    BOOST_CHECK_CLOSE( d, 0.0, 0.01 );
  }
  {
    T const d = box( V::make(  0.0, -1.0, -1.0) );
    BOOST_CHECK_CLOSE( d, 0.0, 0.01 );
  }

  // outside points
  {
    T const d = box( V::make(  2.0,  0.0,  0.0) );
    BOOST_CHECK_CLOSE( d, 1.0, 0.01 );
  }
  {
    T const d = box( V::make(  0.0,  2.0,  0.0) );
    BOOST_CHECK_CLOSE( d, 1.0, 0.01 );
  }
  {
    T const d = box( V::make(  0.0,  0.0,  2.0) );
    BOOST_CHECK_CLOSE( d, 1.0, 0.01 );
  }
  {
    T const d = box( V::make( -2.0,  0.0,  0.0) );
    BOOST_CHECK_CLOSE( d, 1.0, 0.01 );
  }
  {
    T const d = box( V::make(  0.0, -2.0,  0.0) );
    BOOST_CHECK_CLOSE( d, 1.0, 0.01 );
  }
  {
    T const d = box( V::make(  0.0,  0.0, -2.0) );
    BOOST_CHECK_CLOSE( d, 1.0, 0.01 );
  }
  {
    T const d = box( V::make(  2.0,  2.0, 2.0) );
    BOOST_CHECK_CLOSE( d, sqrt(3.0), 0.01 );
  }

  // inside points
  {
    T const d = box( V::make( 0.5,  0.0,  0.0) );
    BOOST_CHECK_CLOSE( d, -0.5, 0.01 );
  }
  {
    T const d = box( V::make(-0.5,  0.0,  0.0) );
    BOOST_CHECK_CLOSE( d, -0.5, 0.01 );
  }
  {
    T const d = box( V::make( 0.0,  0.5,  0.0) );
    BOOST_CHECK_CLOSE( d, -0.5, 0.01 );
  }
  {
    T const d = box( V::make( 0.5, -0.5,  0.0) );
    BOOST_CHECK_CLOSE( d, -0.5, 0.01 );
  }


  //  BOOST_CHECK( ??? );

}

BOOST_AUTO_TEST_SUITE_END();
