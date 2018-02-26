#include <tiny.h>
#include <convex.h>

#define BOOST_AUTO_TEST_MAIN
#include <boost/test/auto_unit_test.hpp>
#include <boost/test/unit_test_suite.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/test_tools.hpp>

BOOST_AUTO_TEST_SUITE(convex_reduce_simplex);

BOOST_AUTO_TEST_CASE(case_by_case_test)
{
  typedef tiny::MathTypes<double>           math_types;
  typedef math_types::vector3_type          vector3_type;
  typedef convex::Simplex<vector3_type>     simplex_type;

  vector3_type const zero = vector3_type::make(0.0, 0.0, 0.0);

  // Simplex is a single point
  {
    vector3_type const v1 = vector3_type::make(  1.0,  1.0,  1.0);

    simplex_type S;

    convex::add_point_to_simplex( v1, v1, v1, S);
    
    vector3_type v;
    vector3_type a;
    vector3_type b;

    v = convex::reduce_simplex( S, a, b );

    BOOST_CHECK( convex::dimension( S ) == 1u );

    BOOST_CHECK( v == zero ); 
    BOOST_CHECK( a == v1 );
    BOOST_CHECK( b == v1 );
  }
  // Simplex is an edge that can not be reduced
  {
    vector3_type const v1 = vector3_type::make(  1.0,  0.0,   1.0);
    vector3_type const v2 = vector3_type::make(  1.0,  0.0,  -1.0);

    simplex_type S;

    convex::add_point_to_simplex( v1, v1, v1, S);
    convex::add_point_to_simplex( v2, v2, v2, S);
    
    vector3_type v;
    vector3_type a;
    vector3_type b;

    v = convex::reduce_simplex( S, a, b );

    BOOST_CHECK( convex::dimension( S ) == 2u );

    BOOST_CHECK( v == zero );

    vector3_type const tst = vector3_type::make(  1.0,  0.0,  0.0);
    BOOST_CHECK( a == tst );
    BOOST_CHECK( b == tst );
  }
  // Simplex is an edge that can be reduced to a single vertex
  {
    vector3_type const v1 = vector3_type::make(  1.0,  1.0,  1.0);
    vector3_type const v2 = vector3_type::make(  2.0,  2.0,  2.0);

    simplex_type S;

    convex::add_point_to_simplex( v1, v1, v1, S);
    convex::add_point_to_simplex( v2, v2, v2, S);
    
    vector3_type v;
    vector3_type a;
    vector3_type b;

    v = convex::reduce_simplex( S, a, b );

    BOOST_CHECK( convex::dimension( S ) == 1u );

    BOOST_CHECK( v == zero );
    BOOST_CHECK( a == v1 );
    BOOST_CHECK( b == v1 );
  }
  // Simplex is a triangle that can not be reduced
  {
    vector3_type const v1 = vector3_type::make(  1.0,  0.0,  0.0);
    vector3_type const v2 = vector3_type::make(  0.0,  1.0,  0.0);
    vector3_type const v3 = vector3_type::make(  0.0,  0.0,  1.0);

    simplex_type S;

    convex::add_point_to_simplex( v1, v1, v1, S);
    convex::add_point_to_simplex( v2, v2, v2, S);
    convex::add_point_to_simplex( v3, v3, v3, S);
    
    vector3_type v;
    vector3_type a;
    vector3_type b;

    v = convex::reduce_simplex( S, a, b );

    BOOST_CHECK( convex::dimension( S ) == 3u );

    BOOST_CHECK( v == zero );
    BOOST_CHECK_CLOSE( a(0), 0.33333333333333331, 0.01);
    BOOST_CHECK_CLOSE( a(1), 0.33333333333333331, 0.01);
    BOOST_CHECK_CLOSE( a(2), 0.33333333333333331, 0.01);
    BOOST_CHECK_CLOSE( b(0), 0.33333333333333331, 0.01);
    BOOST_CHECK_CLOSE( b(1), 0.33333333333333331, 0.01);
    BOOST_CHECK_CLOSE( b(2), 0.33333333333333331, 0.01);    
  }
  // Simplex is a triangle that can be reduced to an edge
  {
    vector3_type const v1 = vector3_type::make(  1.0, -1.0,  0.0);
    vector3_type const v2 = vector3_type::make(  2.0,  0.0,  0.0);
    vector3_type const v3 = vector3_type::make(  1.0,  1.0,  0.0);

    simplex_type S;

    convex::add_point_to_simplex( v1, v1, v1, S);
    convex::add_point_to_simplex( v2, v2, v2, S);
    convex::add_point_to_simplex( v3, v3, v3, S);
    
    vector3_type v;
    vector3_type a;
    vector3_type b;

    v = convex::reduce_simplex( S, a, b );

    BOOST_CHECK( convex::dimension( S ) == 2u );

    BOOST_CHECK( v == zero );
    BOOST_CHECK_CLOSE( a(0), 1.0, 0.01);
    BOOST_CHECK_CLOSE( a(1), 0.0, 0.01);
    BOOST_CHECK_CLOSE( a(2), 0.0, 0.01);
    BOOST_CHECK_CLOSE( b(0), 1.0, 0.01);
    BOOST_CHECK_CLOSE( b(1), 0.0, 0.01);
    BOOST_CHECK_CLOSE( b(2), 0.0, 0.01);
  }
  // Simplex is a triangle that can be reduced to a point
  {
    vector3_type const v1 = vector3_type::make(  2.0, -1.0, 0.0);
    vector3_type const v2 = vector3_type::make(  2.0,  1.0,  0.0);
    vector3_type const v3 = vector3_type::make(  1.0,  0.0,  0.0);

    simplex_type S;

    convex::add_point_to_simplex( v1, v1, v1, S);
    convex::add_point_to_simplex( v2, v2, v2, S);
    convex::add_point_to_simplex( v3, v3, v3, S);
    
    vector3_type v;
    vector3_type a;
    vector3_type b;

    v = convex::reduce_simplex( S, a, b );

    BOOST_CHECK( convex::dimension( S ) == 1u );

    BOOST_CHECK( v == zero );
    BOOST_CHECK_CLOSE( a(0), 1.0, 0.01);
    BOOST_CHECK_CLOSE( a(1), 0.0, 0.01);
    BOOST_CHECK_CLOSE( a(2), 0.0, 0.01);
    BOOST_CHECK_CLOSE( b(0), 1.0, 0.01);
    BOOST_CHECK_CLOSE( b(1), 0.0, 0.01);
    BOOST_CHECK_CLOSE( b(2), 0.0, 0.01);
  }
  // Origin is inside a tetrahedron, so tetrahedron can not be reduced
  {
    vector3_type const v1 = vector3_type::make( -1.0,  0.0, -1.0);
    vector3_type const v2 = vector3_type::make(  1.0,  0.0, -1.0);
    vector3_type const v3 = vector3_type::make(  0.0, -1.0,  1.0);
    vector3_type const v4 = vector3_type::make(  0.0,  1.0,  1.0);

    simplex_type S;

    convex::add_point_to_simplex( v1, v1, v1, S);
    convex::add_point_to_simplex( v2, v2, v2, S);
    convex::add_point_to_simplex( v3, v3, v3, S);
    convex::add_point_to_simplex( v4, v4, v4, S);

    vector3_type v;
    vector3_type a;
    vector3_type b;

    v = convex::reduce_simplex( S, a, b );

    BOOST_CHECK( convex::dimension( S ) == 4u );

    BOOST_CHECK( v == zero );
    BOOST_CHECK( a == v);
    BOOST_CHECK( b == v);
  }
  // Simplex is a tetrahedron that can be reduced to a triangle
  {
    vector3_type const v1 = vector3_type::make(  1.0,  0.0,  0.0);
    vector3_type const v2 = vector3_type::make(  0.0,  1.0,  0.0);
    vector3_type const v3 = vector3_type::make(  0.0,  0.0,  1.0);
    vector3_type const v4 = vector3_type::make(  2.0,  2.0,  2.0);

    simplex_type S;

    convex::add_point_to_simplex( v1, v1, v1, S);
    convex::add_point_to_simplex( v2, v2, v2, S);
    convex::add_point_to_simplex( v3, v3, v3, S);
    convex::add_point_to_simplex( v4, v4, v4, S);

    vector3_type v;
    vector3_type a;
    vector3_type b;

    v = convex::reduce_simplex( S, a, b );

    BOOST_CHECK( convex::dimension( S ) == 3u );

    BOOST_CHECK( v == zero );
    BOOST_CHECK_CLOSE( a(0), 0.33333333333333331, 0.01);
    BOOST_CHECK_CLOSE( a(1), 0.33333333333333331, 0.01);
    BOOST_CHECK_CLOSE( a(2), 0.33333333333333331, 0.01);
    BOOST_CHECK_CLOSE( b(0), 0.33333333333333331, 0.01);
    BOOST_CHECK_CLOSE( b(1), 0.33333333333333331, 0.01);
    BOOST_CHECK_CLOSE( b(2), 0.33333333333333331, 0.01);    
  }
  // Simplex is a tetrahedron that can be reduced to a point
  {
    vector3_type const v1 = vector3_type::make(  1.0,  0.0,  0.0);
    vector3_type const v2 = vector3_type::make(  0.0,  1.0,  0.0);
    vector3_type const v3 = vector3_type::make(  0.0,  0.0,  1.0);
    vector3_type const v4 = vector3_type::make(  0.25,  0.25,  0.25);

    simplex_type S;

    convex::add_point_to_simplex( v1, v1, v1, S);
    convex::add_point_to_simplex( v2, v2, v2, S);
    convex::add_point_to_simplex( v3, v3, v3, S);
    convex::add_point_to_simplex( v4, v4, v4, S);

    vector3_type v;
    vector3_type a;
    vector3_type b;

    v = convex::reduce_simplex( S, a, b );

    BOOST_CHECK( convex::dimension( S ) == 1u );

    BOOST_CHECK( v == zero );
    BOOST_CHECK( a == v4 );
    BOOST_CHECK( b == v4 );
  }
  // Simplex is a tetrahedron that can be reduced to an edge
  {
    vector3_type const v1 = vector3_type::make(  1.0,  -1.0,  0.0);
    vector3_type const v2 = vector3_type::make(  1.0,   1.0,  0.0);
    vector3_type const v3 = vector3_type::make(  2.0,  0.0,  -1.0);
    vector3_type const v4 = vector3_type::make(  2.0,  0.0,   1.0);

    simplex_type S;

    convex::add_point_to_simplex( v1, v1, v1, S);
    convex::add_point_to_simplex( v2, v2, v2, S);
    convex::add_point_to_simplex( v3, v3, v3, S);
    convex::add_point_to_simplex( v4, v4, v4, S);

    vector3_type v;
    vector3_type a;
    vector3_type b;

    v = convex::reduce_simplex( S, a, b );

    BOOST_CHECK( convex::dimension( S ) == 2u );

    BOOST_CHECK( v == zero );
    BOOST_CHECK_CLOSE( a(0), 1.0, 0.01);
    BOOST_CHECK_CLOSE( a(1), 0.0, 0.01);
    BOOST_CHECK_CLOSE( a(2), 0.0, 0.01);
    BOOST_CHECK_CLOSE( b(0), 1.0, 0.01);
    BOOST_CHECK_CLOSE( b(1), 0.0, 0.01);
    BOOST_CHECK_CLOSE( b(2), 0.0, 0.01);
  }
}


BOOST_AUTO_TEST_SUITE_END();
