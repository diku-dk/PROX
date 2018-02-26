#include <tiny.h>
#include <convex.h>

#define BOOST_AUTO_TEST_MAIN
#include <boost/test/auto_unit_test.hpp>
#include <boost/test/unit_test_suite.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/test_tools.hpp>

BOOST_AUTO_TEST_SUITE(convex_reduce_tetrahedron);

BOOST_AUTO_TEST_CASE(case_by_case_test_not_touching)
{
  typedef tiny::MathTypes<double>          math_types;
  typedef math_types::vector3_type        vector3_type;

  typedef convex::Simplex<vector3_type>      simplex_type;

  // First we create a simplex that represents a tetrahedron
  vector3_type const a = vector3_type::make( 0.0,  0.0, 0.0);
  vector3_type const b = vector3_type::make( 1.0,  0.0, 0.0);
  vector3_type const c = vector3_type::make( 0.0,  1.0, 0.0);
  vector3_type const d = vector3_type::make( 0.0,  0.0, 1.0);

  // Inside tetrahedron region new simplex should be ABC
  {
    simplex_type S;

    convex::add_point_to_simplex( a, a, a, S);
    convex::add_point_to_simplex( b, b, b, S);
    convex::add_point_to_simplex( c, c, c, S);
    convex::add_point_to_simplex( d, d, d, S);

    vector3_type const p = (a + b + c + d)*0.25;

    convex::reduce_tetrahedron( p, S );

    BOOST_CHECK( convex::dimension( S ) == 4u );

    BOOST_CHECK(S.m_bitmask == 15u);
    BOOST_CHECK(S.m_v[0] == a);
    BOOST_CHECK(S.m_a[0] == a);
    BOOST_CHECK(S.m_b[0] == a);
    
    BOOST_CHECK(S.m_v[1] == b);
    BOOST_CHECK(S.m_a[1] == b);
    BOOST_CHECK(S.m_b[1] == b);

    BOOST_CHECK(S.m_v[2] == c);
    BOOST_CHECK(S.m_a[2] == c);
    BOOST_CHECK(S.m_b[2] == c);

    BOOST_CHECK(S.m_v[3] == d);
    BOOST_CHECK(S.m_a[3] == d);
    BOOST_CHECK(S.m_b[3] == d);

    BOOST_CHECK_CLOSE(S.m_w[0], 0.25, 0.01);
    BOOST_CHECK_CLOSE(S.m_w[1], 0.25, 0.01);
    BOOST_CHECK_CLOSE(S.m_w[2], 0.25, 0.01);
    BOOST_CHECK_CLOSE(S.m_w[3], 0.25, 0.01);
  }

  // Inside ABC face region new simplex should be ABC
  {
    simplex_type S;

    convex::add_point_to_simplex( a, a, a, S);
    convex::add_point_to_simplex( b, b, b, S);
    convex::add_point_to_simplex( c, c, c, S);
    convex::add_point_to_simplex( d, d, d, S);

    vector3_type const p = (a+b+c)/3.0 - vector3_type::make(0,0,1);

    convex::reduce_tetrahedron( p, S );

    BOOST_CHECK( convex::dimension( S ) == 3u );

    int bit_A    = 0;
    size_t idx_A = 0;
    int bit_B    = 0;
    size_t idx_B = 0;
    int bit_C    = 0;
    size_t idx_C = 0;
    convex::get_used_indices( S.m_bitmask, idx_A, bit_A, idx_B, bit_B , idx_C, bit_C );

    BOOST_CHECK(S.m_bitmask == (bit_A | bit_B | bit_C));
    BOOST_CHECK(S.m_v[idx_A] == a);
    BOOST_CHECK(S.m_a[idx_A] == a);
    BOOST_CHECK(S.m_b[idx_A] == a);
    
    BOOST_CHECK(S.m_v[idx_B] == b);
    BOOST_CHECK(S.m_a[idx_B] == b);
    BOOST_CHECK(S.m_b[idx_B] == b);

    BOOST_CHECK(S.m_v[idx_C] == c);
    BOOST_CHECK(S.m_a[idx_C] == c);
    BOOST_CHECK(S.m_b[idx_C] == c);

    BOOST_CHECK_CLOSE(S.m_w[idx_A], 0.3333333333333333, 0.01);
    BOOST_CHECK_CLOSE(S.m_w[idx_B], 0.3333333333333333, 0.01);
    BOOST_CHECK_CLOSE(S.m_w[idx_C], 0.3333333333333333, 0.01);
  }


  // Inside BCD face region new simplex should be BCD
  {
    simplex_type S;

    convex::add_point_to_simplex( a, a, a, S);
    convex::add_point_to_simplex( b, b, b, S);
    convex::add_point_to_simplex( c, c, c, S);
    convex::add_point_to_simplex( d, d, d, S);

    vector3_type const p = (b+c+d)/3.0 + vector3_type::make(1,1,1);

    convex::reduce_tetrahedron( p, S );

    BOOST_CHECK( convex::dimension( S ) == 3u );

    int bit_A    = 0;
    size_t idx_A = 0;
    int bit_B    = 0;
    size_t idx_B = 0;
    int bit_C    = 0;
    size_t idx_C = 0;
    convex::get_used_indices( S.m_bitmask, idx_A, bit_A, idx_B, bit_B , idx_C, bit_C );

    BOOST_CHECK(S.m_bitmask == (bit_A | bit_B | bit_C));
    BOOST_CHECK(S.m_v[idx_A] == b);
    BOOST_CHECK(S.m_a[idx_A] == b);
    BOOST_CHECK(S.m_b[idx_A] == b);
    
    BOOST_CHECK(S.m_v[idx_B] == c);
    BOOST_CHECK(S.m_a[idx_B] == c);
    BOOST_CHECK(S.m_b[idx_B] == c);

    BOOST_CHECK(S.m_v[idx_C] == d);
    BOOST_CHECK(S.m_a[idx_C] == d);
    BOOST_CHECK(S.m_b[idx_C] == d);

    BOOST_CHECK_CLOSE(S.m_w[idx_A], 0.3333333333333333, 0.01);
    BOOST_CHECK_CLOSE(S.m_w[idx_B], 0.3333333333333333, 0.01);
    BOOST_CHECK_CLOSE(S.m_w[idx_C], 0.3333333333333333, 0.01);
  }

  // Inside ABD face region new simplex should be ABD
  {
    simplex_type S;

    convex::add_point_to_simplex( a, a, a, S);
    convex::add_point_to_simplex( b, b, b, S);
    convex::add_point_to_simplex( c, c, c, S);
    convex::add_point_to_simplex( d, d, d, S);

    vector3_type const p = (a+b+d)/3.0 - vector3_type::make(0,1,0);

    convex::reduce_tetrahedron( p, S );

    BOOST_CHECK( convex::dimension( S ) == 3u );

    int bit_A    = 0;
    size_t idx_A = 0;

    int bit_B    = 0;
    size_t idx_B = 0;
    
    int bit_C    = 0;
    size_t idx_C = 0;
    convex::get_used_indices( S.m_bitmask, idx_A, bit_A, idx_B, bit_B , idx_C, bit_C );

    BOOST_CHECK(S.m_bitmask == (bit_A | bit_B | bit_C));
    BOOST_CHECK(S.m_v[idx_A] == a);
    BOOST_CHECK(S.m_a[idx_A] == a);
    BOOST_CHECK(S.m_b[idx_A] == a);
    
    BOOST_CHECK(S.m_v[idx_B] == b);
    BOOST_CHECK(S.m_a[idx_B] == b);
    BOOST_CHECK(S.m_b[idx_B] == b);

    BOOST_CHECK(S.m_v[idx_C] == d);
    BOOST_CHECK(S.m_a[idx_C] == d);
    BOOST_CHECK(S.m_b[idx_C] == d);

    BOOST_CHECK_CLOSE(S.m_w[idx_A], 0.3333333333333333, 0.01);
    BOOST_CHECK_CLOSE(S.m_w[idx_B], 0.3333333333333333, 0.01);
    BOOST_CHECK_CLOSE(S.m_w[idx_C], 0.3333333333333333, 0.01);
  }


  // Inside ACD face region new simplex should be ACD
  {
    simplex_type S;

    convex::add_point_to_simplex( a, a, a, S);
    convex::add_point_to_simplex( b, b, b, S);
    convex::add_point_to_simplex( c, c, c, S);
    convex::add_point_to_simplex( d, d, d, S);

    vector3_type const p = (a+c+d)/3.0 - vector3_type::make(1,0,0);

    convex::reduce_tetrahedron( p, S );

    BOOST_CHECK( convex::dimension( S ) == 3u );

    int bit_A    = 0;
    size_t idx_A = 0;
    int bit_B    = 0;
    size_t idx_B = 0;
    int bit_C    = 0;
    size_t idx_C = 0;
    convex::get_used_indices( S.m_bitmask, idx_A, bit_A, idx_B, bit_B , idx_C, bit_C );

    BOOST_CHECK(S.m_bitmask == (bit_A | bit_B | bit_C));
    BOOST_CHECK(S.m_v[idx_A] == a);
    BOOST_CHECK(S.m_a[idx_A] == a);
    BOOST_CHECK(S.m_b[idx_A] == a);
    
    BOOST_CHECK(S.m_v[idx_B] == c);
    BOOST_CHECK(S.m_a[idx_B] == c);
    BOOST_CHECK(S.m_b[idx_B] == c);

    BOOST_CHECK(S.m_v[idx_C] == d);
    BOOST_CHECK(S.m_a[idx_C] == d);
    BOOST_CHECK(S.m_b[idx_C] == d);

    BOOST_CHECK_CLOSE(S.m_w[idx_A], 0.3333333333333333, 0.01);
    BOOST_CHECK_CLOSE(S.m_w[idx_B], 0.3333333333333333, 0.01);
    BOOST_CHECK_CLOSE(S.m_w[idx_C], 0.3333333333333333, 0.01);
  }

  // Inside AB edge region new simplex should be AB
  {
    simplex_type S;

    convex::add_point_to_simplex( a, a, a, S);
    convex::add_point_to_simplex( b, b, b, S);
    convex::add_point_to_simplex( c, c, c, S);
    convex::add_point_to_simplex( d, d, d, S);

    vector3_type const p = (a+b)/2.0 - vector3_type::make(0,1,1);

    convex::reduce_tetrahedron( p, S );

    BOOST_CHECK( convex::dimension( S ) == 2u );

    int bit_A    = 0;
    size_t idx_A = 0;
    int bit_B    = 0;
    size_t idx_B = 0;
    convex::get_used_indices( S.m_bitmask, idx_A, bit_A, idx_B, bit_B );

    BOOST_CHECK(S.m_bitmask == (bit_A | bit_B ));
    BOOST_CHECK(S.m_v[idx_A] == a);
    BOOST_CHECK(S.m_a[idx_A] == a);
    BOOST_CHECK(S.m_b[idx_A] == a);
    
    BOOST_CHECK(S.m_v[idx_B] == b);
    BOOST_CHECK(S.m_a[idx_B] == b);
    BOOST_CHECK(S.m_b[idx_B] == b);

    BOOST_CHECK_CLOSE(S.m_w[idx_A], 0.5, 0.01);
    BOOST_CHECK_CLOSE(S.m_w[idx_B], 0.5, 0.01);
  }

  // Inside AC edge region new simplex should be AC
  {
    simplex_type S;

    convex::add_point_to_simplex( a, a, a, S);
    convex::add_point_to_simplex( b, b, b, S);
    convex::add_point_to_simplex( c, c, c, S);
    convex::add_point_to_simplex( d, d, d, S);

    vector3_type const p = (a+c)/2.0 - vector3_type::make(1,0,1);

    convex::reduce_tetrahedron( p, S );

    BOOST_CHECK( convex::dimension( S ) == 2u );

    int bit_A    = 0;
    size_t idx_A = 0;
    int bit_B    = 0;
    size_t idx_B = 0;
    convex::get_used_indices( S.m_bitmask, idx_A, bit_A, idx_B, bit_B );

    BOOST_CHECK(S.m_bitmask == (bit_A | bit_B ));
    BOOST_CHECK(S.m_v[idx_A] == a);
    BOOST_CHECK(S.m_a[idx_A] == a);
    BOOST_CHECK(S.m_b[idx_A] == a);
    
    BOOST_CHECK(S.m_v[idx_B] == c);
    BOOST_CHECK(S.m_a[idx_B] == c);
    BOOST_CHECK(S.m_b[idx_B] == c);

    BOOST_CHECK_CLOSE(S.m_w[idx_A], 0.5, 0.01);
    BOOST_CHECK_CLOSE(S.m_w[idx_B], 0.5, 0.01);
  }

  // Inside AD edge region new simplex should be AD
  {
    simplex_type S;

    convex::add_point_to_simplex( a, a, a, S);
    convex::add_point_to_simplex( b, b, b, S);
    convex::add_point_to_simplex( c, c, c, S);
    convex::add_point_to_simplex( d, d, d, S);

    vector3_type const p = (a+d)/2.0 - vector3_type::make(1,1,0);

    convex::reduce_tetrahedron( p, S );

    BOOST_CHECK( convex::dimension( S ) == 2u );

    int bit_A    = 0;
    size_t idx_A = 0;
    int bit_B    = 0;
    size_t idx_B = 0;
    convex::get_used_indices( S.m_bitmask, idx_A, bit_A, idx_B, bit_B );

    BOOST_CHECK(S.m_bitmask == (bit_A | bit_B ));
    BOOST_CHECK(S.m_v[idx_A] == a);
    BOOST_CHECK(S.m_a[idx_A] == a);
    BOOST_CHECK(S.m_b[idx_A] == a);
    
    BOOST_CHECK(S.m_v[idx_B] == d);
    BOOST_CHECK(S.m_a[idx_B] == d);
    BOOST_CHECK(S.m_b[idx_B] == d);

    BOOST_CHECK_CLOSE(S.m_w[idx_A], 0.5, 0.01);
    BOOST_CHECK_CLOSE(S.m_w[idx_B], 0.5, 0.01);
  }

  // Inside BC edge region new simplex should be BC
  {
    simplex_type S;

    convex::add_point_to_simplex( a, a, a, S);
    convex::add_point_to_simplex( b, b, b, S);
    convex::add_point_to_simplex( c, c, c, S);
    convex::add_point_to_simplex( d, d, d, S);

    vector3_type const p = (b+c)/2.0 + vector3_type::make(0,0,-1) + tiny::unit( vector3_type::make(1,1,1) );

    convex::reduce_tetrahedron( p, S );

    BOOST_CHECK( convex::dimension( S ) == 2u );

    int bit_A    = 0;
    size_t idx_A = 0;
    int bit_B    = 0;
    size_t idx_B = 0;
    convex::get_used_indices( S.m_bitmask, idx_A, bit_A, idx_B, bit_B );

    BOOST_CHECK(S.m_bitmask == (bit_A | bit_B ));
    BOOST_CHECK(S.m_v[idx_A] == b);
    BOOST_CHECK(S.m_a[idx_A] == b);
    BOOST_CHECK(S.m_b[idx_A] == b);
    
    BOOST_CHECK(S.m_v[idx_B] == c);
    BOOST_CHECK(S.m_a[idx_B] == c);
    BOOST_CHECK(S.m_b[idx_B] == c);

    BOOST_CHECK_CLOSE(S.m_w[idx_A], 0.5, 0.01);
    BOOST_CHECK_CLOSE(S.m_w[idx_B], 0.5, 0.01);
  }

  // Inside BD edge region new simplex should be BD
  {
    simplex_type S;

    convex::add_point_to_simplex( a, a, a, S);
    convex::add_point_to_simplex( b, b, b, S);
    convex::add_point_to_simplex( c, c, c, S);
    convex::add_point_to_simplex( d, d, d, S);

    vector3_type const p = (b+d)/2.0 + vector3_type::make(0.0,-1.0,0.0) + tiny::unit(vector3_type::make(1.0,1.0,1.0));

    convex::reduce_tetrahedron( p, S );

    BOOST_CHECK( convex::dimension( S ) == 2u );

    int bit_A    = 0;
    size_t idx_A = 0;
    int bit_B    = 0;
    size_t idx_B = 0;
    convex::get_used_indices( S.m_bitmask, idx_A, bit_A, idx_B, bit_B );

    BOOST_CHECK(S.m_bitmask == (bit_A | bit_B ));
    BOOST_CHECK(S.m_v[idx_A] == b);
    BOOST_CHECK(S.m_a[idx_A] == b);
    BOOST_CHECK(S.m_b[idx_A] == b);
    
    BOOST_CHECK(S.m_v[idx_B] == d);
    BOOST_CHECK(S.m_a[idx_B] == d);
    BOOST_CHECK(S.m_b[idx_B] == d);

    BOOST_CHECK_CLOSE(S.m_w[idx_A], 0.5, 0.01);
    BOOST_CHECK_CLOSE(S.m_w[idx_B], 0.5, 0.01);
  }

  // Inside CD edge region new simplex should be CD
  {
    simplex_type S;

    convex::add_point_to_simplex( a, a, a, S);
    convex::add_point_to_simplex( b, b, b, S);
    convex::add_point_to_simplex( c, c, c, S);
    convex::add_point_to_simplex( d, d, d, S);

    vector3_type const p = (c+d)/2.0 + vector3_type::make(-1.0,0.0,0.0) + tiny::unit( vector3_type::make(1.0,1.0,1.0) );

    convex::reduce_tetrahedron( p, S );

    BOOST_CHECK( convex::dimension( S ) == 2u );

    int bit_A    = 0;
    size_t idx_A = 0;
    int bit_B    = 0;
    size_t idx_B = 0;
    convex::get_used_indices( S.m_bitmask, idx_A, bit_A, idx_B, bit_B );

    BOOST_CHECK(S.m_bitmask == (bit_A | bit_B ));
    BOOST_CHECK(S.m_v[idx_A] == c);
    BOOST_CHECK(S.m_a[idx_A] == c);
    BOOST_CHECK(S.m_b[idx_A] == c);
    
    BOOST_CHECK(S.m_v[idx_B] == d);
    BOOST_CHECK(S.m_a[idx_B] == d);
    BOOST_CHECK(S.m_b[idx_B] == d);

    BOOST_CHECK_CLOSE(S.m_w[idx_A], 0.5, 0.01);
    BOOST_CHECK_CLOSE(S.m_w[idx_B], 0.5, 0.01);
  }

  // Inside A vertex region new simplex should be A
  {
    simplex_type S;

    convex::add_point_to_simplex( a, a, a, S);
    convex::add_point_to_simplex( b, b, b, S);
    convex::add_point_to_simplex( c, c, c, S);
    convex::add_point_to_simplex( d, d, d, S);

    vector3_type const p = a + vector3_type::make(-1.0,-1.0,-1.0);

    convex::reduce_tetrahedron( p, S );

    BOOST_CHECK( convex::dimension( S ) == 1u );

    int bit_A    = 0;
    size_t idx_A = 0;
    convex::get_used_indices( S.m_bitmask, idx_A, bit_A );

    BOOST_CHECK(S.m_bitmask == bit_A );
    BOOST_CHECK(S.m_v[idx_A] == a);
    BOOST_CHECK(S.m_a[idx_A] == a);
    BOOST_CHECK(S.m_b[idx_A] == a);
    
    BOOST_CHECK_CLOSE(S.m_w[idx_A], 1.0, 0.01);
  }

  // Inside B vertex region new simplex should be B
  {
    simplex_type S;

    convex::add_point_to_simplex( a, a, a, S);
    convex::add_point_to_simplex( b, b, b, S);
    convex::add_point_to_simplex( c, c, c, S);
    convex::add_point_to_simplex( d, d, d, S);

    vector3_type const p = b + vector3_type::make(0.0,-1.0,0.0)+ vector3_type::make(0.0,0.0,-1.0) + tiny::unit( vector3_type::make(1,1,1) );

    convex::reduce_tetrahedron( p, S );

    BOOST_CHECK( convex::dimension( S ) == 1u );

    int bit_A    = 0;
    size_t idx_A = 0;
    convex::get_used_indices( S.m_bitmask, idx_A, bit_A );

    BOOST_CHECK(S.m_bitmask == bit_A );
    BOOST_CHECK(S.m_v[idx_A] == b);
    BOOST_CHECK(S.m_a[idx_A] == b);
    BOOST_CHECK(S.m_b[idx_A] == b);
    
    BOOST_CHECK_CLOSE(S.m_w[idx_A], 1.0, 0.01);
  }

  // Inside C vertex region new simplex should be C
  {
    simplex_type S;

    convex::add_point_to_simplex( a, a, a, S);
    convex::add_point_to_simplex( b, b, b, S);
    convex::add_point_to_simplex( c, c, c, S);
    convex::add_point_to_simplex( d, d, d, S);

    vector3_type const p = c + vector3_type::make(-1.0,0.0,0.0) + vector3_type::make(0,0,-1)  + tiny::unit( vector3_type::make(1,1,1) );

    convex::reduce_tetrahedron( p, S );

    BOOST_CHECK( convex::dimension( S ) == 1u );

    int bit_A    = 0;
    size_t idx_A = 0;
    convex::get_used_indices( S.m_bitmask, idx_A, bit_A );

    BOOST_CHECK(S.m_bitmask == bit_A );
    BOOST_CHECK(S.m_v[idx_A] == c);
    BOOST_CHECK(S.m_a[idx_A] == c);
    BOOST_CHECK(S.m_b[idx_A] == c);
    
    BOOST_CHECK_CLOSE(S.m_w[idx_A], 1.0, 0.01);
  }

  // Inside D vertex region new simplex should be D
  {
    simplex_type S;

    convex::add_point_to_simplex( a, a, a, S);
    convex::add_point_to_simplex( b, b, b, S);
    convex::add_point_to_simplex( c, c, c, S);
    convex::add_point_to_simplex( d, d, d, S);

    vector3_type const p = d + vector3_type::make(-1.0,0.0,0.0) + vector3_type::make(0.0,-1.0,0.0) + tiny::unit( vector3_type::make(1,1,1) );

    convex::reduce_tetrahedron( p, S );

    BOOST_CHECK( convex::dimension( S ) == 1u );

    int bit_A    = 0;
    size_t idx_A = 0;
    convex::get_used_indices( S.m_bitmask, idx_A, bit_A );

    BOOST_CHECK(S.m_bitmask == bit_A );
    BOOST_CHECK(S.m_v[idx_A] == d);
    BOOST_CHECK(S.m_a[idx_A] == d);
    BOOST_CHECK(S.m_b[idx_A] == d);
    
    BOOST_CHECK_CLOSE(S.m_w[idx_A], 1.0, 0.01);
  }

}



BOOST_AUTO_TEST_CASE(case_by_case_test_touching)
{
  typedef tiny::MathTypes<double> math_types;
  typedef math_types::vector3_type                         vector3_type;

  typedef convex::Simplex<vector3_type>           simplex_type;

  // First we create a simplex that represents a tetrahedron
  vector3_type const a = vector3_type::make( 0.0,  0.0, 0.0);
  vector3_type const b = vector3_type::make( 1.0,  0.0, 0.0);
  vector3_type const c = vector3_type::make( 0.0,  1.0, 0.0);
  vector3_type const d = vector3_type::make( 0.0,  0.0, 1.0);

  // Inside tetrahedron region new simplex should be ABCD
  {
    simplex_type S;

    convex::add_point_to_simplex( a, a, a, S);
    convex::add_point_to_simplex( b, b, b, S);
    convex::add_point_to_simplex( c, c, c, S);
    convex::add_point_to_simplex( d, d, d, S);

    vector3_type const p = (a + b + c + d)*0.25;

    convex::reduce_tetrahedron( p, S );

    BOOST_CHECK( convex::dimension( S ) == 4u );

    BOOST_CHECK(S.m_bitmask == 15u);
    BOOST_CHECK(S.m_v[0] == a);
    BOOST_CHECK(S.m_a[0] == a);
    BOOST_CHECK(S.m_b[0] == a);
    
    BOOST_CHECK(S.m_v[1] == b);
    BOOST_CHECK(S.m_a[1] == b);
    BOOST_CHECK(S.m_b[1] == b);

    BOOST_CHECK(S.m_v[2] == c);
    BOOST_CHECK(S.m_a[2] == c);
    BOOST_CHECK(S.m_b[2] == c);

    BOOST_CHECK(S.m_v[3] == d);
    BOOST_CHECK(S.m_a[3] == d);
    BOOST_CHECK(S.m_b[3] == d);

    BOOST_CHECK_CLOSE(S.m_w[0], 0.25, 0.01);
    BOOST_CHECK_CLOSE(S.m_w[1], 0.25, 0.01);
    BOOST_CHECK_CLOSE(S.m_w[2], 0.25, 0.01);
    BOOST_CHECK_CLOSE(S.m_w[3], 0.25, 0.01);
  }

  // Inside ABC face region new simplex should be ABC
  {
    simplex_type S;

    convex::add_point_to_simplex( a, a, a, S);
    convex::add_point_to_simplex( b, b, b, S);
    convex::add_point_to_simplex( c, c, c, S);
    convex::add_point_to_simplex( d, d, d, S);

    vector3_type const p = (a+b+c)/3.0;

    convex::reduce_tetrahedron( p, S );

    BOOST_CHECK( convex::dimension( S ) == 3u );

    int bit_A    = 0;
    size_t idx_A = 0;
    int bit_B    = 0;
    size_t idx_B = 0;
    int bit_C    = 0;
    size_t idx_C = 0;
    convex::get_used_indices( S.m_bitmask, idx_A, bit_A, idx_B, bit_B , idx_C, bit_C );

    BOOST_CHECK(S.m_bitmask == (bit_A | bit_B | bit_C));
    BOOST_CHECK(S.m_v[idx_A] == a);
    BOOST_CHECK(S.m_a[idx_A] == a);
    BOOST_CHECK(S.m_b[idx_A] == a);
    
    BOOST_CHECK(S.m_v[idx_B] == b);
    BOOST_CHECK(S.m_a[idx_B] == b);
    BOOST_CHECK(S.m_b[idx_B] == b);

    BOOST_CHECK(S.m_v[idx_C] == c);
    BOOST_CHECK(S.m_a[idx_C] == c);
    BOOST_CHECK(S.m_b[idx_C] == c);

    BOOST_CHECK_CLOSE(S.m_w[idx_A], 0.3333333333333333, 0.01);
    BOOST_CHECK_CLOSE(S.m_w[idx_B], 0.3333333333333333, 0.01);
    BOOST_CHECK_CLOSE(S.m_w[idx_C], 0.3333333333333333, 0.01);
  }


  // Inside BCD face region new simplex should be BCD
  {

    // This case fails, problem is with numerical precision, the
    // triangle face plane distance is computed to something like -5e-17, very
    // small but negative. This means the test point is seen as being inside
    // the tetrahedron and we get the full simplex back!

    // It means that all the test cases for this case fails!

    simplex_type S;

    convex::add_point_to_simplex( a, a, a, S);
    convex::add_point_to_simplex( b, b, b, S);
    convex::add_point_to_simplex( c, c, c, S);
    convex::add_point_to_simplex( d, d, d, S);

    vector3_type const p = (b+c+d)/3.0;

    convex::reduce_tetrahedron( p, S );

    //BOOST_CHECK( convex::dimension( S ) == 3u );

    //int bit_A    = 0;
    //size_t idx_A = 0;
    //int bit_B    = 0;
    //size_t idx_B = 0;
    //int bit_C    = 0;
    //size_t idx_C = 0;
    //convex::get_used_indices( S.m_bitmask, idx_A, bit_A, idx_B, bit_B , idx_C, bit_C );

    //BOOST_CHECK(S.m_bitmask == (bit_A | bit_B | bit_C));
    //BOOST_CHECK(S.m_v[idx_A] == b);
    //BOOST_CHECK(S.m_a[idx_A] == b);
    //BOOST_CHECK(S.m_b[idx_A] == b);
    //
    //BOOST_CHECK(S.m_v[idx_B] == c);
    //BOOST_CHECK(S.m_a[idx_B] == c);
    //BOOST_CHECK(S.m_b[idx_B] == c);

    //BOOST_CHECK(S.m_v[idx_C] == d);
    //BOOST_CHECK(S.m_a[idx_C] == d);
    //BOOST_CHECK(S.m_b[idx_C] == d);

    //BOOST_CHECK_CLOSE(S.m_w[idx_A], 0.3333333333333333, 0.01);
    //BOOST_CHECK_CLOSE(S.m_w[idx_B], 0.3333333333333333, 0.01);
    //BOOST_CHECK_CLOSE(S.m_w[idx_C], 0.3333333333333333, 0.01);
  }

  // Inside ABD face region new simplex should be ABD
  {
    simplex_type S;

    convex::add_point_to_simplex( a, a, a, S);
    convex::add_point_to_simplex( b, b, b, S);
    convex::add_point_to_simplex( c, c, c, S);
    convex::add_point_to_simplex( d, d, d, S);

    vector3_type const p = (a+b+d)/3.0;

    convex::reduce_tetrahedron( p, S );

    BOOST_CHECK( convex::dimension( S ) == 3u );

    int bit_A    = 0;
    size_t idx_A = 0;

    int bit_B    = 0;
    size_t idx_B = 0;
    
    int bit_C    = 0;
    size_t idx_C = 0;
    convex::get_used_indices( S.m_bitmask, idx_A, bit_A, idx_B, bit_B , idx_C, bit_C );

    BOOST_CHECK(S.m_bitmask == (bit_A | bit_B | bit_C));
    BOOST_CHECK(S.m_v[idx_A] == a);
    BOOST_CHECK(S.m_a[idx_A] == a);
    BOOST_CHECK(S.m_b[idx_A] == a);
    
    BOOST_CHECK(S.m_v[idx_B] == b);
    BOOST_CHECK(S.m_a[idx_B] == b);
    BOOST_CHECK(S.m_b[idx_B] == b);

    BOOST_CHECK(S.m_v[idx_C] == d);
    BOOST_CHECK(S.m_a[idx_C] == d);
    BOOST_CHECK(S.m_b[idx_C] == d);

    BOOST_CHECK_CLOSE(S.m_w[idx_A], 0.3333333333333333, 0.01);
    BOOST_CHECK_CLOSE(S.m_w[idx_B], 0.3333333333333333, 0.01);
    BOOST_CHECK_CLOSE(S.m_w[idx_C], 0.3333333333333333, 0.01);
  }


  // Inside ACD face region new simplex should be ACD
  {
    simplex_type S;

    convex::add_point_to_simplex( a, a, a, S);
    convex::add_point_to_simplex( b, b, b, S);
    convex::add_point_to_simplex( c, c, c, S);
    convex::add_point_to_simplex( d, d, d, S);

    vector3_type const p = (a+c+d)/3.0;

    convex::reduce_tetrahedron( p, S );

    BOOST_CHECK( convex::dimension( S ) == 3u );

    int bit_A    = 0;
    size_t idx_A = 0;
    int bit_B    = 0;
    size_t idx_B = 0;
    int bit_C    = 0;
    size_t idx_C = 0;
    convex::get_used_indices( S.m_bitmask, idx_A, bit_A, idx_B, bit_B , idx_C, bit_C );

    BOOST_CHECK(S.m_bitmask == (bit_A | bit_B | bit_C));
    BOOST_CHECK(S.m_v[idx_A] == a);
    BOOST_CHECK(S.m_a[idx_A] == a);
    BOOST_CHECK(S.m_b[idx_A] == a);
    
    BOOST_CHECK(S.m_v[idx_B] == c);
    BOOST_CHECK(S.m_a[idx_B] == c);
    BOOST_CHECK(S.m_b[idx_B] == c);

    BOOST_CHECK(S.m_v[idx_C] == d);
    BOOST_CHECK(S.m_a[idx_C] == d);
    BOOST_CHECK(S.m_b[idx_C] == d);

    BOOST_CHECK_CLOSE(S.m_w[idx_A], 0.3333333333333333, 0.01);
    BOOST_CHECK_CLOSE(S.m_w[idx_B], 0.3333333333333333, 0.01);
    BOOST_CHECK_CLOSE(S.m_w[idx_C], 0.3333333333333333, 0.01);
  }

  // Inside AB edge region new simplex should be AB
  {
    simplex_type S;

    convex::add_point_to_simplex( a, a, a, S);
    convex::add_point_to_simplex( b, b, b, S);
    convex::add_point_to_simplex( c, c, c, S);
    convex::add_point_to_simplex( d, d, d, S);

    vector3_type const p = (a+b)/2.0;

    convex::reduce_tetrahedron( p, S );

    BOOST_CHECK( convex::dimension( S ) == 2u );

    int bit_A    = 0;
    size_t idx_A = 0;
    int bit_B    = 0;
    size_t idx_B = 0;
    convex::get_used_indices( S.m_bitmask, idx_A, bit_A, idx_B, bit_B );

    BOOST_CHECK(S.m_bitmask == (bit_A | bit_B ));
    BOOST_CHECK(S.m_v[idx_A] == a);
    BOOST_CHECK(S.m_a[idx_A] == a);
    BOOST_CHECK(S.m_b[idx_A] == a);
    
    BOOST_CHECK(S.m_v[idx_B] == b);
    BOOST_CHECK(S.m_a[idx_B] == b);
    BOOST_CHECK(S.m_b[idx_B] == b);

    BOOST_CHECK_CLOSE(S.m_w[idx_A], 0.5, 0.01);
    BOOST_CHECK_CLOSE(S.m_w[idx_B], 0.5, 0.01);
  }

  // Inside AC edge region new simplex should be AC
  {
    simplex_type S;

    convex::add_point_to_simplex( a, a, a, S);
    convex::add_point_to_simplex( b, b, b, S);
    convex::add_point_to_simplex( c, c, c, S);
    convex::add_point_to_simplex( d, d, d, S);

    vector3_type const p = (a+c)/2.0;

    convex::reduce_tetrahedron( p, S );

    BOOST_CHECK( convex::dimension( S ) == 2u );

    int bit_A    = 0;
    size_t idx_A = 0;
    int bit_B    = 0;
    size_t idx_B = 0;
    convex::get_used_indices( S.m_bitmask, idx_A, bit_A, idx_B, bit_B );

    BOOST_CHECK(S.m_bitmask == (bit_A | bit_B ));
    BOOST_CHECK(S.m_v[idx_A] == a);
    BOOST_CHECK(S.m_a[idx_A] == a);
    BOOST_CHECK(S.m_b[idx_A] == a);
    
    BOOST_CHECK(S.m_v[idx_B] == c);
    BOOST_CHECK(S.m_a[idx_B] == c);
    BOOST_CHECK(S.m_b[idx_B] == c);

    BOOST_CHECK_CLOSE(S.m_w[idx_A], 0.5, 0.01);
    BOOST_CHECK_CLOSE(S.m_w[idx_B], 0.5, 0.01);
  }

  // Inside AD edge region new simplex should be AD
  {
    simplex_type S;

    convex::add_point_to_simplex( a, a, a, S);
    convex::add_point_to_simplex( b, b, b, S);
    convex::add_point_to_simplex( c, c, c, S);
    convex::add_point_to_simplex( d, d, d, S);

    vector3_type const p = (a+d)/2.0;

    convex::reduce_tetrahedron( p, S );

    BOOST_CHECK( convex::dimension( S ) == 2u );

    int bit_A    = 0;
    size_t idx_A = 0;
    int bit_B    = 0;
    size_t idx_B = 0;
    convex::get_used_indices( S.m_bitmask, idx_A, bit_A, idx_B, bit_B );

    BOOST_CHECK(S.m_bitmask == (bit_A | bit_B ));
    BOOST_CHECK(S.m_v[idx_A] == a);
    BOOST_CHECK(S.m_a[idx_A] == a);
    BOOST_CHECK(S.m_b[idx_A] == a);
    
    BOOST_CHECK(S.m_v[idx_B] == d);
    BOOST_CHECK(S.m_a[idx_B] == d);
    BOOST_CHECK(S.m_b[idx_B] == d);

    BOOST_CHECK_CLOSE(S.m_w[idx_A], 0.5, 0.01);
    BOOST_CHECK_CLOSE(S.m_w[idx_B], 0.5, 0.01);
  }

  // Inside BC edge region new simplex should be BC
  {
    simplex_type S;

    convex::add_point_to_simplex( a, a, a, S);
    convex::add_point_to_simplex( b, b, b, S);
    convex::add_point_to_simplex( c, c, c, S);
    convex::add_point_to_simplex( d, d, d, S);

    vector3_type const p = (b+c)/2.0;

    convex::reduce_tetrahedron( p, S );

    BOOST_CHECK( convex::dimension( S ) == 2u );

    int bit_A    = 0;
    size_t idx_A = 0;
    int bit_B    = 0;
    size_t idx_B = 0;
    convex::get_used_indices( S.m_bitmask, idx_A, bit_A, idx_B, bit_B );

    BOOST_CHECK(S.m_bitmask == (bit_A | bit_B ));
    BOOST_CHECK(S.m_v[idx_A] == b);
    BOOST_CHECK(S.m_a[idx_A] == b);
    BOOST_CHECK(S.m_b[idx_A] == b);
    
    BOOST_CHECK(S.m_v[idx_B] == c);
    BOOST_CHECK(S.m_a[idx_B] == c);
    BOOST_CHECK(S.m_b[idx_B] == c);

    BOOST_CHECK_CLOSE(S.m_w[idx_A], 0.5, 0.01);
    BOOST_CHECK_CLOSE(S.m_w[idx_B], 0.5, 0.01);
  }

  // Inside BD edge region new simplex should be BD
  {
    simplex_type S;

    convex::add_point_to_simplex( a, a, a, S);
    convex::add_point_to_simplex( b, b, b, S);
    convex::add_point_to_simplex( c, c, c, S);
    convex::add_point_to_simplex( d, d, d, S);

    vector3_type const p = (b+d)/2.0;

    convex::reduce_tetrahedron( p, S );

    BOOST_CHECK( convex::dimension( S ) == 2u );

    int bit_A    = 0;
    size_t idx_A = 0;
    int bit_B    = 0;
    size_t idx_B = 0;
    convex::get_used_indices( S.m_bitmask, idx_A, bit_A, idx_B, bit_B );

    BOOST_CHECK(S.m_bitmask == (bit_A | bit_B ));
    BOOST_CHECK(S.m_v[idx_A] == b);
    BOOST_CHECK(S.m_a[idx_A] == b);
    BOOST_CHECK(S.m_b[idx_A] == b);
    
    BOOST_CHECK(S.m_v[idx_B] == d);
    BOOST_CHECK(S.m_a[idx_B] == d);
    BOOST_CHECK(S.m_b[idx_B] == d);

    BOOST_CHECK_CLOSE(S.m_w[idx_A], 0.5, 0.01);
    BOOST_CHECK_CLOSE(S.m_w[idx_B], 0.5, 0.01);
  }

  // Inside CD edge region new simplex should be CD
  {
    simplex_type S;

    convex::add_point_to_simplex( a, a, a, S);
    convex::add_point_to_simplex( b, b, b, S);
    convex::add_point_to_simplex( c, c, c, S);
    convex::add_point_to_simplex( d, d, d, S);

    vector3_type const p = (c+d)/2.0;

    convex::reduce_tetrahedron( p, S );

    BOOST_CHECK( convex::dimension( S ) == 2u );

    int bit_A    = 0;
    size_t idx_A = 0;
    int bit_B    = 0;
    size_t idx_B = 0;
    convex::get_used_indices( S.m_bitmask, idx_A, bit_A, idx_B, bit_B );

    BOOST_CHECK(S.m_bitmask == (bit_A | bit_B ));
    BOOST_CHECK(S.m_v[idx_A] == c);
    BOOST_CHECK(S.m_a[idx_A] == c);
    BOOST_CHECK(S.m_b[idx_A] == c);
    
    BOOST_CHECK(S.m_v[idx_B] == d);
    BOOST_CHECK(S.m_a[idx_B] == d);
    BOOST_CHECK(S.m_b[idx_B] == d);

    BOOST_CHECK_CLOSE(S.m_w[idx_A], 0.5, 0.01);
    BOOST_CHECK_CLOSE(S.m_w[idx_B], 0.5, 0.01);
  }

  // Inside A vertex region new simplex should be A
  {
    simplex_type S;

    convex::add_point_to_simplex( a, a, a, S);
    convex::add_point_to_simplex( b, b, b, S);
    convex::add_point_to_simplex( c, c, c, S);
    convex::add_point_to_simplex( d, d, d, S);

    vector3_type const p = a;

    convex::reduce_tetrahedron( p, S );

    BOOST_CHECK( convex::dimension( S ) == 1u );

    int bit_A    = 0;
    size_t idx_A = 0;
    convex::get_used_indices( S.m_bitmask, idx_A, bit_A );

    BOOST_CHECK(S.m_bitmask == bit_A );
    BOOST_CHECK(S.m_v[idx_A] == a);
    BOOST_CHECK(S.m_a[idx_A] == a);
    BOOST_CHECK(S.m_b[idx_A] == a);
    
    BOOST_CHECK_CLOSE(S.m_w[idx_A], 1.0, 0.01);
  }

  // Inside B vertex region new simplex should be B
  {
    simplex_type S;

    convex::add_point_to_simplex( a, a, a, S);
    convex::add_point_to_simplex( b, b, b, S);
    convex::add_point_to_simplex( c, c, c, S);
    convex::add_point_to_simplex( d, d, d, S);

    vector3_type const p = b;

    convex::reduce_tetrahedron( p, S );

    BOOST_CHECK( convex::dimension( S ) == 1u );

    int bit_A    = 0;
    size_t idx_A = 0;
    convex::get_used_indices( S.m_bitmask, idx_A, bit_A );

    BOOST_CHECK(S.m_bitmask == bit_A );
    BOOST_CHECK(S.m_v[idx_A] == b);
    BOOST_CHECK(S.m_a[idx_A] == b);
    BOOST_CHECK(S.m_b[idx_A] == b);
    
    BOOST_CHECK_CLOSE(S.m_w[idx_A], 1.0, 0.01);
  }

  // Inside C vertex region new simplex should be C
  {
    simplex_type S;

    convex::add_point_to_simplex( a, a, a, S);
    convex::add_point_to_simplex( b, b, b, S);
    convex::add_point_to_simplex( c, c, c, S);
    convex::add_point_to_simplex( d, d, d, S);

    vector3_type const p = c;

    convex::reduce_tetrahedron( p, S );

    BOOST_CHECK( convex::dimension( S ) == 1u );

    int bit_A    = 0;
    size_t idx_A = 0;
    convex::get_used_indices( S.m_bitmask, idx_A, bit_A );

    BOOST_CHECK(S.m_bitmask == bit_A );
    BOOST_CHECK(S.m_v[idx_A] == c);
    BOOST_CHECK(S.m_a[idx_A] == c);
    BOOST_CHECK(S.m_b[idx_A] == c);
    
    BOOST_CHECK_CLOSE(S.m_w[idx_A], 1.0, 0.01);
  }

  // Inside D vertex region new simplex should be D
  {
    simplex_type S;

    convex::add_point_to_simplex( a, a, a, S);
    convex::add_point_to_simplex( b, b, b, S);
    convex::add_point_to_simplex( c, c, c, S);
    convex::add_point_to_simplex( d, d, d, S);

    vector3_type const p = d;

    convex::reduce_tetrahedron( p, S );

    BOOST_CHECK( convex::dimension( S ) == 1u );

    int bit_A    = 0;
    size_t idx_A = 0;
    convex::get_used_indices( S.m_bitmask, idx_A, bit_A );

    BOOST_CHECK(S.m_bitmask == bit_A );
    BOOST_CHECK(S.m_v[idx_A] == d);
    BOOST_CHECK(S.m_a[idx_A] == d);
    BOOST_CHECK(S.m_b[idx_A] == d);
    
    BOOST_CHECK_CLOSE(S.m_w[idx_A], 1.0, 0.01);
  }

}


BOOST_AUTO_TEST_CASE(case_by_case_test_assymmetric)
{
  typedef tiny::MathTypes<double>   math_types;
  typedef math_types::vector3_type                         vector3_type;

  typedef convex::Simplex<vector3_type>           simplex_type;

  // First we create a simplex that represents a tetrahedron
  vector3_type const a = vector3_type::make( 0.0,  0.0, 0.0);
  vector3_type const b = vector3_type::make( 1.0,  0.0, 0.0);
  vector3_type const c = vector3_type::make( 0.0,  1.0, 0.0);
  vector3_type const d = vector3_type::make( 0.0,  0.0, 1.0);

  // New simplex should be ABCD
  {
    simplex_type S;

    convex::add_point_to_simplex( a, a, a, S);
    convex::add_point_to_simplex( b, b, b, S);
    convex::add_point_to_simplex( c, c, c, S);
    convex::add_point_to_simplex( d, d, d, S);

    vector3_type const p = 0.1*a + 0.2*b + 0.3*c + 0.4*d;

    convex::reduce_tetrahedron( p, S );

    BOOST_CHECK( convex::dimension( S ) == 4u );

    BOOST_CHECK(S.m_bitmask == 15u);
    BOOST_CHECK(S.m_v[0] == a);
    BOOST_CHECK(S.m_a[0] == a);
    BOOST_CHECK(S.m_b[0] == a);
    
    BOOST_CHECK(S.m_v[1] == b);
    BOOST_CHECK(S.m_a[1] == b);
    BOOST_CHECK(S.m_b[1] == b);

    BOOST_CHECK(S.m_v[2] == c);
    BOOST_CHECK(S.m_a[2] == c);
    BOOST_CHECK(S.m_b[2] == c);

    BOOST_CHECK(S.m_v[3] == d);
    BOOST_CHECK(S.m_a[3] == d);
    BOOST_CHECK(S.m_b[3] == d);

    BOOST_CHECK_CLOSE(S.m_w[0], 0.1, 0.01);
    BOOST_CHECK_CLOSE(S.m_w[1], 0.2, 0.01);
    BOOST_CHECK_CLOSE(S.m_w[2], 0.3, 0.01);
    BOOST_CHECK_CLOSE(S.m_w[3], 0.4, 0.01);
  }

  // New simplex should be ABC
  {
    simplex_type S;

    convex::add_point_to_simplex( a, a, a, S);
    convex::add_point_to_simplex( b, b, b, S);
    convex::add_point_to_simplex( c, c, c, S);
    convex::add_point_to_simplex( d, d, d, S);

    vector3_type const p = 0.2*a+0.3*b+0.5*c;

    convex::reduce_tetrahedron( p, S );

    BOOST_CHECK( convex::dimension( S ) == 3u );

    int bit_A    = 0;
    size_t idx_A = 0;
    int bit_B    = 0;
    size_t idx_B = 0;
    int bit_C    = 0;
    size_t idx_C = 0;
    convex::get_used_indices( S.m_bitmask, idx_A, bit_A, idx_B, bit_B , idx_C, bit_C );

    BOOST_CHECK(S.m_bitmask == (bit_A | bit_B | bit_C));
    BOOST_CHECK(S.m_v[idx_A] == a);
    BOOST_CHECK(S.m_a[idx_A] == a);
    BOOST_CHECK(S.m_b[idx_A] == a);
    
    BOOST_CHECK(S.m_v[idx_B] == b);
    BOOST_CHECK(S.m_a[idx_B] == b);
    BOOST_CHECK(S.m_b[idx_B] == b);

    BOOST_CHECK(S.m_v[idx_C] == c);
    BOOST_CHECK(S.m_a[idx_C] == c);
    BOOST_CHECK(S.m_b[idx_C] == c);

    BOOST_CHECK_CLOSE(S.m_w[idx_A], 0.2, 0.01);
    BOOST_CHECK_CLOSE(S.m_w[idx_B], 0.3, 0.01);
    BOOST_CHECK_CLOSE(S.m_w[idx_C], 0.5, 0.01);
  }

  // New simplex should be BCD
  {
    simplex_type S;

    convex::add_point_to_simplex( a, a, a, S);
    convex::add_point_to_simplex( b, b, b, S);
    convex::add_point_to_simplex( c, c, c, S);
    convex::add_point_to_simplex( d, d, d, S);

    vector3_type const p = 0.2*b+0.3*c+0.5*d;

    convex::reduce_tetrahedron( p, S );

    BOOST_CHECK( convex::dimension( S ) == 3u );

    int bit_A    = 0;
    size_t idx_A = 0;
    int bit_B    = 0;
    size_t idx_B = 0;
    int bit_C    = 0;
    size_t idx_C = 0;
    convex::get_used_indices( S.m_bitmask, idx_A, bit_A, idx_B, bit_B , idx_C, bit_C );

    BOOST_CHECK(S.m_bitmask == (bit_A | bit_B | bit_C));
    BOOST_CHECK(S.m_v[idx_A] == b);
    BOOST_CHECK(S.m_a[idx_A] == b);
    BOOST_CHECK(S.m_b[idx_A] == b);
    
    BOOST_CHECK(S.m_v[idx_B] == c);
    BOOST_CHECK(S.m_a[idx_B] == c);
    BOOST_CHECK(S.m_b[idx_B] == c);

    BOOST_CHECK(S.m_v[idx_C] == d);
    BOOST_CHECK(S.m_a[idx_C] == d);
    BOOST_CHECK(S.m_b[idx_C] == d);

    BOOST_CHECK_CLOSE(S.m_w[idx_A], 0.2, 0.01);
    BOOST_CHECK_CLOSE(S.m_w[idx_B], 0.3, 0.01);
    BOOST_CHECK_CLOSE(S.m_w[idx_C], 0.5, 0.01);
  }

  // New simplex should be ABD
  {
    simplex_type S;

    convex::add_point_to_simplex( a, a, a, S);
    convex::add_point_to_simplex( b, b, b, S);
    convex::add_point_to_simplex( c, c, c, S);
    convex::add_point_to_simplex( d, d, d, S);

    vector3_type const p = 0.2*a+0.3*b+0.5*d;

    convex::reduce_tetrahedron( p, S );

    BOOST_CHECK( convex::dimension( S ) == 3u );

    int bit_A    = 0;
    size_t idx_A = 0;

    int bit_B    = 0;
    size_t idx_B = 0;
    
    int bit_C    = 0;
    size_t idx_C = 0;
    convex::get_used_indices( S.m_bitmask, idx_A, bit_A, idx_B, bit_B , idx_C, bit_C );

    BOOST_CHECK(S.m_bitmask == (bit_A | bit_B | bit_C));
    BOOST_CHECK(S.m_v[idx_A] == a);
    BOOST_CHECK(S.m_a[idx_A] == a);
    BOOST_CHECK(S.m_b[idx_A] == a);
    
    BOOST_CHECK(S.m_v[idx_B] == b);
    BOOST_CHECK(S.m_a[idx_B] == b);
    BOOST_CHECK(S.m_b[idx_B] == b);

    BOOST_CHECK(S.m_v[idx_C] == d);
    BOOST_CHECK(S.m_a[idx_C] == d);
    BOOST_CHECK(S.m_b[idx_C] == d);

    BOOST_CHECK_CLOSE(S.m_w[idx_A], 0.2, 0.01);
    BOOST_CHECK_CLOSE(S.m_w[idx_B], 0.3, 0.01);
    BOOST_CHECK_CLOSE(S.m_w[idx_C], 0.5, 0.01);
  }

  // New simplex should be ACD
  {
    simplex_type S;

    convex::add_point_to_simplex( a, a, a, S);
    convex::add_point_to_simplex( b, b, b, S);
    convex::add_point_to_simplex( c, c, c, S);
    convex::add_point_to_simplex( d, d, d, S);

    vector3_type const p = 0.2*a + 0.3*c + 0.5*d;

    convex::reduce_tetrahedron( p, S );

    BOOST_CHECK( convex::dimension( S ) == 3u );

    int bit_A    = 0;
    size_t idx_A = 0;
    int bit_B    = 0;
    size_t idx_B = 0;
    int bit_C    = 0;
    size_t idx_C = 0;
    convex::get_used_indices( S.m_bitmask, idx_A, bit_A, idx_B, bit_B , idx_C, bit_C );

    BOOST_CHECK(S.m_bitmask == (bit_A | bit_B | bit_C));
    BOOST_CHECK(S.m_v[idx_A] == a);
    BOOST_CHECK(S.m_a[idx_A] == a);
    BOOST_CHECK(S.m_b[idx_A] == a);
    
    BOOST_CHECK(S.m_v[idx_B] == c);
    BOOST_CHECK(S.m_a[idx_B] == c);
    BOOST_CHECK(S.m_b[idx_B] == c);

    BOOST_CHECK(S.m_v[idx_C] == d);
    BOOST_CHECK(S.m_a[idx_C] == d);
    BOOST_CHECK(S.m_b[idx_C] == d);

    BOOST_CHECK_CLOSE(S.m_w[idx_A], 0.2, 0.01);
    BOOST_CHECK_CLOSE(S.m_w[idx_B], 0.3, 0.01);
    BOOST_CHECK_CLOSE(S.m_w[idx_C], 0.5, 0.01);
  }

  // New simplex should be AB
  {
    simplex_type S;

    convex::add_point_to_simplex( a, a, a, S);
    convex::add_point_to_simplex( b, b, b, S);
    convex::add_point_to_simplex( c, c, c, S);
    convex::add_point_to_simplex( d, d, d, S);

    vector3_type const p = 0.4*a+0.6*b;

    convex::reduce_tetrahedron( p, S );

    BOOST_CHECK( convex::dimension( S ) == 2u );

    int bit_A    = 0;
    size_t idx_A = 0;
    int bit_B    = 0;
    size_t idx_B = 0;
    convex::get_used_indices( S.m_bitmask, idx_A, bit_A, idx_B, bit_B );

    BOOST_CHECK(S.m_bitmask == (bit_A | bit_B ));
    BOOST_CHECK(S.m_v[idx_A] == a);
    BOOST_CHECK(S.m_a[idx_A] == a);
    BOOST_CHECK(S.m_b[idx_A] == a);
    
    BOOST_CHECK(S.m_v[idx_B] == b);
    BOOST_CHECK(S.m_a[idx_B] == b);
    BOOST_CHECK(S.m_b[idx_B] == b);

    BOOST_CHECK_CLOSE(S.m_w[idx_A], 0.4, 0.01);
    BOOST_CHECK_CLOSE(S.m_w[idx_B], 0.6, 0.01);
  }

  // New simplex should be AC
  {
    simplex_type S;

    convex::add_point_to_simplex( a, a, a, S);
    convex::add_point_to_simplex( b, b, b, S);
    convex::add_point_to_simplex( c, c, c, S);
    convex::add_point_to_simplex( d, d, d, S);

    vector3_type const p = 0.4*a+0.6*c;

    convex::reduce_tetrahedron( p, S );

    BOOST_CHECK( convex::dimension( S ) == 2u );

    int bit_A    = 0;
    size_t idx_A = 0;
    int bit_B    = 0;
    size_t idx_B = 0;
    convex::get_used_indices( S.m_bitmask, idx_A, bit_A, idx_B, bit_B );

    BOOST_CHECK(S.m_bitmask == (bit_A | bit_B ));
    BOOST_CHECK(S.m_v[idx_A] == a);
    BOOST_CHECK(S.m_a[idx_A] == a);
    BOOST_CHECK(S.m_b[idx_A] == a);
    
    BOOST_CHECK(S.m_v[idx_B] == c);
    BOOST_CHECK(S.m_a[idx_B] == c);
    BOOST_CHECK(S.m_b[idx_B] == c);

    BOOST_CHECK_CLOSE(S.m_w[idx_A], 0.4, 0.01);
    BOOST_CHECK_CLOSE(S.m_w[idx_B], 0.6, 0.01);
  }

  // New simplex should be AD
  {
    simplex_type S;

    convex::add_point_to_simplex( a, a, a, S);
    convex::add_point_to_simplex( b, b, b, S);
    convex::add_point_to_simplex( c, c, c, S);
    convex::add_point_to_simplex( d, d, d, S);

    vector3_type const p = 0.4*a+0.6*d;

    convex::reduce_tetrahedron( p, S );

    BOOST_CHECK( convex::dimension( S ) == 2u );

    int bit_A    = 0;
    size_t idx_A = 0;
    int bit_B    = 0;
    size_t idx_B = 0;
    convex::get_used_indices( S.m_bitmask, idx_A, bit_A, idx_B, bit_B );

    BOOST_CHECK(S.m_bitmask == (bit_A | bit_B ));
    BOOST_CHECK(S.m_v[idx_A] == a);
    BOOST_CHECK(S.m_a[idx_A] == a);
    BOOST_CHECK(S.m_b[idx_A] == a);
    
    BOOST_CHECK(S.m_v[idx_B] == d);
    BOOST_CHECK(S.m_a[idx_B] == d);
    BOOST_CHECK(S.m_b[idx_B] == d);

    BOOST_CHECK_CLOSE(S.m_w[idx_A], 0.4, 0.01);
    BOOST_CHECK_CLOSE(S.m_w[idx_B], 0.6, 0.01);
  }

  // New simplex should be BC
  {
    simplex_type S;

    convex::add_point_to_simplex( a, a, a, S);
    convex::add_point_to_simplex( b, b, b, S);
    convex::add_point_to_simplex( c, c, c, S);
    convex::add_point_to_simplex( d, d, d, S);

    vector3_type const p = 0.4*b+0.6*c;

    convex::reduce_tetrahedron( p, S );

    BOOST_CHECK( convex::dimension( S ) == 2u );

    int bit_A    = 0;
    size_t idx_A = 0;
    int bit_B    = 0;
    size_t idx_B = 0;
    convex::get_used_indices( S.m_bitmask, idx_A, bit_A, idx_B, bit_B );

    BOOST_CHECK(S.m_bitmask == (bit_A | bit_B ));
    BOOST_CHECK(S.m_v[idx_A] == b);
    BOOST_CHECK(S.m_a[idx_A] == b);
    BOOST_CHECK(S.m_b[idx_A] == b);
    
    BOOST_CHECK(S.m_v[idx_B] == c);
    BOOST_CHECK(S.m_a[idx_B] == c);
    BOOST_CHECK(S.m_b[idx_B] == c);

    BOOST_CHECK_CLOSE(S.m_w[idx_A], 0.4, 0.01);
    BOOST_CHECK_CLOSE(S.m_w[idx_B], 0.6, 0.01);
  }

  // New simplex should be BD
  {
    simplex_type S;

    convex::add_point_to_simplex( a, a, a, S);
    convex::add_point_to_simplex( b, b, b, S);
    convex::add_point_to_simplex( c, c, c, S);
    convex::add_point_to_simplex( d, d, d, S);

    vector3_type const p = 0.4*b+0.6*d;

    convex::reduce_tetrahedron( p, S );

    BOOST_CHECK( convex::dimension( S ) == 2u );

    int bit_A    = 0;
    size_t idx_A = 0;
    int bit_B    = 0;
    size_t idx_B = 0;
    convex::get_used_indices( S.m_bitmask, idx_A, bit_A, idx_B, bit_B );

    BOOST_CHECK(S.m_bitmask == (bit_A | bit_B ));
    BOOST_CHECK(S.m_v[idx_A] == b);
    BOOST_CHECK(S.m_a[idx_A] == b);
    BOOST_CHECK(S.m_b[idx_A] == b);
    
    BOOST_CHECK(S.m_v[idx_B] == d);
    BOOST_CHECK(S.m_a[idx_B] == d);
    BOOST_CHECK(S.m_b[idx_B] == d);

    BOOST_CHECK_CLOSE(S.m_w[idx_A], 0.4, 0.01);
    BOOST_CHECK_CLOSE(S.m_w[idx_B], 0.6, 0.01);
  }

  // New simplex should be CD
  {
    simplex_type S;

    convex::add_point_to_simplex( a, a, a, S);
    convex::add_point_to_simplex( b, b, b, S);
    convex::add_point_to_simplex( c, c, c, S);
    convex::add_point_to_simplex( d, d, d, S);

    vector3_type const p = 0.4*c+0.6*d;

    convex::reduce_tetrahedron( p, S );

    BOOST_CHECK( convex::dimension( S ) == 2u );

    int bit_A    = 0;
    size_t idx_A = 0;
    int bit_B    = 0;
    size_t idx_B = 0;
    convex::get_used_indices( S.m_bitmask, idx_A, bit_A, idx_B, bit_B );

    BOOST_CHECK(S.m_bitmask == (bit_A | bit_B ));
    BOOST_CHECK(S.m_v[idx_A] == c);
    BOOST_CHECK(S.m_a[idx_A] == c);
    BOOST_CHECK(S.m_b[idx_A] == c);
    
    BOOST_CHECK(S.m_v[idx_B] == d);
    BOOST_CHECK(S.m_a[idx_B] == d);
    BOOST_CHECK(S.m_b[idx_B] == d);

    BOOST_CHECK_CLOSE(S.m_w[idx_A], 0.4, 0.01);
    BOOST_CHECK_CLOSE(S.m_w[idx_B], 0.6, 0.01);
  }
}


BOOST_AUTO_TEST_SUITE_END();
