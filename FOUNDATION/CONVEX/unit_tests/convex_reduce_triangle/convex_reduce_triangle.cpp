#include <tiny.h>
#include <convex.h>

#define BOOST_AUTO_TEST_MAIN
#include <boost/test/auto_unit_test.hpp>
#include <boost/test/unit_test_suite.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/test_tools.hpp>


// 2011-11-12 Kenny: Does these test depend on the order in which vertices are added to the simplex? We always use a,b,c

BOOST_AUTO_TEST_SUITE(convex_reduce_triangle);

BOOST_AUTO_TEST_CASE(case_by_case_test)
{
  typedef tiny::MathTypes<double>         math_types;
  typedef math_types::vector3_type        vector3_type;
  typedef math_types::real_type           real_type;
  
  typedef convex::Simplex<vector3_type>           simplex_type;
  
  // Inside face-region new simplex should be ABC
  {
    vector3_type const a = vector3_type::make(-1.0, -1.0, 0.0);
    vector3_type const b = vector3_type::make( 1.0, -1.0, 0.0);
    vector3_type const c = vector3_type::make( 0.5,  0.1, 0.0);
    
    simplex_type S;
    
    convex::add_point_to_simplex( a, a, a, S);
    convex::add_point_to_simplex( b, b, b, S);
    convex::add_point_to_simplex( c, c, c, S);
    
    for( size_t i=1u; i<10u; ++i)
      for( size_t j=1u; j<10u; ++j)
        for( size_t k=1u; k<10u; ++k)
        {
          real_type u = i*0.1;
          real_type v = j*0.1;
          real_type w = k*0.1;
          
          real_type lgh = u + v + w;
          u /= lgh;
          v /= lgh;
          w /= lgh;
          
          vector3_type const p = u*a + v*b + w*c;
          
          convex::reduce_triangle( p, S );
          
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
          BOOST_CHECK_CLOSE(S.m_w[idx_A], u, 0.01);
          BOOST_CHECK_CLOSE(S.m_w[idx_B], v, 0.01);
          BOOST_CHECK_CLOSE(S.m_w[idx_C], w, 0.01);
          
        }
  }
  // First we create a simplex that represents an triangle
  vector3_type const a = vector3_type::make(-1.0, -1.0, 0.0);
  vector3_type const b = vector3_type::make( 1.0, -1.0, 0.0);
  vector3_type const c = vector3_type::make( 0.0,  1.0, 0.0);
  
  // Inside face-region new simplex should be ABC
  {
    simplex_type S;
    
    convex::add_point_to_simplex( a, a, a, S);
    convex::add_point_to_simplex( b, b, b, S);
    convex::add_point_to_simplex( c, c, c, S);
    
    vector3_type const p = vector3_type::make( 0.0, 0.0,  0.0);
    
    convex::reduce_triangle( p, S );
    
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
    
    BOOST_CHECK_CLOSE(S.m_w[idx_A], 0.25, 0.01);
    BOOST_CHECK_CLOSE(S.m_w[idx_B], 0.25, 0.01);
    BOOST_CHECK_CLOSE(S.m_w[idx_C], 0.5, 0.01);
  }
  // Inside A voronoi region new simplex should be A
  {
    simplex_type S;
    
    convex::add_point_to_simplex( a, a, a, S);
    convex::add_point_to_simplex( b, b, b, S);
    convex::add_point_to_simplex( c, c, c, S);
    
    vector3_type const p = vector3_type::make( -2.0, -1.0,  1.0);
    
    convex::reduce_triangle( p, S );
    
    BOOST_CHECK( convex::dimension( S ) == 1u );
    
    int bit_A    = 0;
    size_t idx_A = 0;
    convex::get_used_indices( S.m_bitmask, idx_A, bit_A );
    
    BOOST_CHECK(S.m_bitmask == bit_A);
    BOOST_CHECK(S.m_v[idx_A] == a);
    BOOST_CHECK(S.m_a[idx_A] == a);
    BOOST_CHECK(S.m_b[idx_A] == a);
    
    BOOST_CHECK_CLOSE(S.m_w[idx_A], 1.0, 0.01);
  }
  // Inside B voronoi region new simplex should be B
  {
    simplex_type S;
    
    convex::add_point_to_simplex( a, a, a, S);
    convex::add_point_to_simplex( b, b, b, S);
    convex::add_point_to_simplex( c, c, c, S);
    
    vector3_type const p = vector3_type::make(  2.0, -1.0,  1.0);
    
    convex::reduce_triangle( p, S );
    
    BOOST_CHECK( convex::dimension( S ) == 1u );
    
    int bit_A    = 0;
    size_t idx_A = 0;
    convex::get_used_indices( S.m_bitmask, idx_A, bit_A );
    
    BOOST_CHECK(S.m_bitmask == bit_A);
    BOOST_CHECK(S.m_v[idx_A] == b);
    BOOST_CHECK(S.m_a[idx_A] == b);
    BOOST_CHECK(S.m_b[idx_A] == b);
    
    BOOST_CHECK_CLOSE(S.m_w[idx_A], 1.00, 0.01);
  }
  // Inside C voronoi region new simplex should be C
  {
    simplex_type S;
    
    convex::add_point_to_simplex( a, a, a, S);
    convex::add_point_to_simplex( b, b, b, S);
    convex::add_point_to_simplex( c, c, c, S);
    
    vector3_type const p = vector3_type::make(  0.0, 2.0,  1.0);
    
    convex::reduce_triangle( p, S );
    
    BOOST_CHECK( convex::dimension( S ) == 1u );
    
    int bit_A    = 0;
    size_t idx_A = 0;
    convex::get_used_indices( S.m_bitmask, idx_A, bit_A  );
    
    BOOST_CHECK(S.m_bitmask == bit_A);
    BOOST_CHECK(S.m_v[idx_A] == c);
    BOOST_CHECK(S.m_a[idx_A] == c);
    BOOST_CHECK(S.m_b[idx_A] == c);
    
    BOOST_CHECK_CLOSE(S.m_w[idx_A], 1.0, 0.01);
  }
  // Inside AB voronoi region new simplex should be AB
  {
    simplex_type S;
    
    convex::add_point_to_simplex( a, a, a, S);
    convex::add_point_to_simplex( b, b, b, S);
    convex::add_point_to_simplex( c, c, c, S);
    
    vector3_type const p = vector3_type::make( 0.0, -2.0,  1.0);
    
    convex::reduce_triangle( p, S );
    
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
  // Inside BC voronoi region new simplex should be BC
  {
    simplex_type S;
    
    convex::add_point_to_simplex( a, a, a, S);
    convex::add_point_to_simplex( b, b, b, S);
    convex::add_point_to_simplex( c, c, c, S);
    
    vector3_type const p = vector3_type::make( 1.5, 0.5,  1.0);
    
    convex::reduce_triangle( p, S );
    
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
  // Inside AC voronoi region new simplex should be AC
  {
    simplex_type S;
    
    convex::add_point_to_simplex( a, a, a, S);
    convex::add_point_to_simplex( b, b, b, S);
    convex::add_point_to_simplex( c, c, c, S);
    
    vector3_type const p = vector3_type::make( -1.5, 0.5,  1.0);
    
    convex::reduce_triangle( p, S );
    
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
  // On vertex A new simplex should be A
  {
    simplex_type S;
    
    convex::add_point_to_simplex( a, a, a, S);
    convex::add_point_to_simplex( b, b, b, S);
    convex::add_point_to_simplex( c, c, c, S);
    
    convex::reduce_triangle( a, S );
    
    BOOST_CHECK( convex::dimension( S ) == 1u );
    
    int bit_A    = 0;
    size_t idx_A = 0;
    convex::get_used_indices( S.m_bitmask, idx_A, bit_A );
    
    BOOST_CHECK(S.m_bitmask == bit_A);
    BOOST_CHECK(S.m_v[idx_A] == a);
    BOOST_CHECK(S.m_a[idx_A] == a);
    BOOST_CHECK(S.m_b[idx_A] == a);
    
    BOOST_CHECK_CLOSE(S.m_w[idx_A], 1.0, 0.01);
  }
  // On vertex B new simplex should be B
  {
    simplex_type S;
    
    convex::add_point_to_simplex( a, a, a, S);
    convex::add_point_to_simplex( b, b, b, S);
    convex::add_point_to_simplex( c, c, c, S);
    
    convex::reduce_triangle( b, S );
    
    BOOST_CHECK( convex::dimension( S ) == 1u );
    
    int bit_A    = 0;
    size_t idx_A = 0;
    convex::get_used_indices( S.m_bitmask, idx_A, bit_A );
    
    BOOST_CHECK(S.m_bitmask == bit_A);
    BOOST_CHECK(S.m_v[idx_A] == b);
    BOOST_CHECK(S.m_a[idx_A] == b);
    BOOST_CHECK(S.m_b[idx_A] == b);
    
    BOOST_CHECK_CLOSE(S.m_w[idx_A], 1.00, 0.01);
  }
  // On vertex C new simplex should be C
  {
    simplex_type S;
    
    convex::add_point_to_simplex( a, a, a, S);
    convex::add_point_to_simplex( b, b, b, S);
    convex::add_point_to_simplex( c, c, c, S);
    
    convex::reduce_triangle( c, S );
    
    BOOST_CHECK( convex::dimension( S ) == 1u );
    
    int bit_A    = 0;
    size_t idx_A = 0;
    convex::get_used_indices( S.m_bitmask, idx_A, bit_A  );
    
    BOOST_CHECK(S.m_bitmask == bit_A);
    BOOST_CHECK(S.m_v[idx_A] == c);
    BOOST_CHECK(S.m_a[idx_A] == c);
    BOOST_CHECK(S.m_b[idx_A] == c);
    
    BOOST_CHECK_CLOSE(S.m_w[idx_A], 1.0, 0.01);
  }
  // On edge AB new simplex should be AB
  {
    simplex_type S;
    
    convex::add_point_to_simplex( a, a, a, S);
    convex::add_point_to_simplex( b, b, b, S);
    convex::add_point_to_simplex( c, c, c, S);
    
    vector3_type const p = vector3_type::make( 0.0, -1.0,  1.0);
    
    convex::reduce_triangle( p, S );
    
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
  // On edge BC new simplex should be BC
  {
    simplex_type S;
    
    convex::add_point_to_simplex( a, a, a, S);
    convex::add_point_to_simplex( b, b, b, S);
    convex::add_point_to_simplex( c, c, c, S);
    
    vector3_type const p = vector3_type::make( 0.5, 0.0,  1.0);// this point is inside the triangle?!
    
    convex::reduce_triangle( p, S );
    
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
  // On edge AC new simplex should be AC
  {
    simplex_type S;
    
    convex::add_point_to_simplex( a, a, a, S);
    convex::add_point_to_simplex( b, b, b, S);
    convex::add_point_to_simplex( c, c, c, S);
    
    vector3_type const p = vector3_type::make( -0.5, 0.0,  1.0);
    
    convex::reduce_triangle( p, S );
    
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
  // Assymmetric test cases
  
  // New simplex should be AB
  {
    simplex_type S;
    
    convex::add_point_to_simplex( a, a, a, S);
    convex::add_point_to_simplex( b, b, b, S);
    convex::add_point_to_simplex( c, c, c, S);
    
    vector3_type const p = 0.4*a + 0.6*b;
    
    convex::reduce_triangle( p, S );
    
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
  // New simplex should be BC
  {
    simplex_type S;
    
    convex::add_point_to_simplex( a, a, a, S);
    convex::add_point_to_simplex( b, b, b, S);
    convex::add_point_to_simplex( c, c, c, S);
    
    vector3_type const p = 0.4*b + 0.6*c;// point is inside triangle?!
    
    convex::reduce_triangle( p, S );
    
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
  // New simplex should be AC
  {
    simplex_type S;
    
    convex::add_point_to_simplex( a, a, a, S);
    convex::add_point_to_simplex( b, b, b, S);
    convex::add_point_to_simplex( c, c, c, S);
    
    vector3_type const p = 0.6*a + 0.4*c;
    
    convex::reduce_triangle( p, S );
    
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
    
    BOOST_CHECK_CLOSE(S.m_w[idx_A], 0.6, 0.01);
    BOOST_CHECK_CLOSE(S.m_w[idx_B], 0.4, 0.01);
  }
  // New simplex should be ABC
  {
    simplex_type S;
    
    convex::add_point_to_simplex( a, a, a, S);
    convex::add_point_to_simplex( b, b, b, S);
    convex::add_point_to_simplex( c, c, c, S);
    
    vector3_type const p = 0.1*a + 0.2*b + 0.7*c;
    
    convex::reduce_triangle( p, S );
    
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
    
    BOOST_CHECK_CLOSE(S.m_w[idx_A], 0.1, 0.01);
    BOOST_CHECK_CLOSE(S.m_w[idx_B], 0.2, 0.01);
    BOOST_CHECK_CLOSE(S.m_w[idx_C], 0.7, 0.01);
  }
}

BOOST_AUTO_TEST_SUITE_END();
