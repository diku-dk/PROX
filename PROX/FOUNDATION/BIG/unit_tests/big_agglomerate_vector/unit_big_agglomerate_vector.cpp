#include <big_types.h>
#include <min_map/big_agglomerate_vector.h>

#define BOOST_AUTO_TEST_MAIN
#include <boost/test/auto_unit_test.hpp>
#include <boost/test/unit_test_suite.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/test_tools.hpp>

BOOST_AUTO_TEST_SUITE(big_agglomerate_vector);

BOOST_AUTO_TEST_CASE(test_case1)
{
  typedef ublas::vector<size_t>       idx_vector_type;
  typedef ublas::vector<double>       vector_type;

  idx_vector_type new2old;
  idx_vector_type old2new;

  old2new.resize(10,false);
  new2old.resize(10,false);

  old2new( 2 ) = 0 ;
  old2new( 4 ) = 1 ;
  old2new( 6 ) = 2 ;
  old2new( 9 ) = 3 ;
  old2new( 1 ) = 4 ;
  old2new( 3 ) = 5 ;
  old2new( 8 ) = 6 ;
  old2new( 0 ) = 7 ;
  old2new( 5 ) = 8 ;
  old2new( 7 ) = 9 ;

  new2old( 0 ) = 2 ;
  new2old( 1 ) = 4 ;
  new2old( 2 ) = 6 ;
  new2old( 3 ) = 9 ;
  new2old( 4 ) = 1 ;
  new2old( 5 ) = 3 ;
  new2old( 6 ) = 8 ;
  new2old( 7 ) = 0 ;
  new2old( 8 ) = 5 ;
  new2old( 9 ) = 7 ;

  vector_type x_a;
  vector_type x_b;
  vector_type x;
  x.resize(10,false);

  x_a.resize(4,false);
  x_b.resize(6,false);

  x_a(0) = 1.0;
  x_a(1) = 2.0;
  x_a(2) = 3.0;
  x_a(3) = 4.0;

  x_b(0) = 5.0;
  x_b(1) = 6.0;
  x_b(2) = 7.0;
  x_b(3) = 8.0;
  x_b(4) = 9.0;
  x_b(5) = 10.0;

  x.clear();

  big::agglomerate_vector( x_a, x_b, new2old, x);

  double tol = 0.01;

  BOOST_CHECK_CLOSE( double( x( 2 ) ), double( x_a(0) ), tol );
  BOOST_CHECK_CLOSE( double( x( 4 ) ), double( x_a(1) ), tol );
  BOOST_CHECK_CLOSE( double( x( 6 ) ), double( x_a(2) ), tol );
  BOOST_CHECK_CLOSE( double( x( 9 ) ), double( x_a(3) ), tol );
  BOOST_CHECK_CLOSE( double( x( 1 ) ), double( x_b(0) ), tol );
  BOOST_CHECK_CLOSE( double( x( 3 ) ), double( x_b(1) ), tol );
  BOOST_CHECK_CLOSE( double( x( 8 ) ), double( x_b(2) ), tol );
  BOOST_CHECK_CLOSE( double( x( 0 ) ), double( x_b(3) ), tol );
  BOOST_CHECK_CLOSE( double( x( 5 ) ), double( x_b(4) ), tol );
  BOOST_CHECK_CLOSE( double( x( 7 ) ), double( x_b(5) ), tol );
}

BOOST_AUTO_TEST_CASE(test_case2)
{
  typedef ublas::vector<size_t>       idx_vector_type;
  typedef ublas::vector<double>       vector_type;
  
  idx_vector_type new2old;
  idx_vector_type old2new;
  
  old2new.resize(10,false);
  new2old.resize(10,false);
  
  old2new( 0 ) = 0 ;
  old2new( 1 ) = 1 ;
  old2new( 2 ) = 2 ;
  old2new( 3 ) = 3 ;
  old2new( 4 ) = 4 ;
  old2new( 5 ) = 5 ;
  old2new( 6 ) = 6 ;
  old2new( 7 ) = 7 ;
  old2new( 8 ) = 8 ;
  old2new( 9 ) = 9 ;
  
  new2old( 0 ) = 0 ;
  new2old( 1 ) = 1 ;
  new2old( 2 ) = 2 ;
  new2old( 3 ) = 3 ;
  new2old( 4 ) = 4 ;
  new2old( 5 ) = 5 ;
  new2old( 6 ) = 6 ;
  new2old( 7 ) = 7 ;
  new2old( 8 ) = 8 ;
  new2old( 9 ) = 9 ;
  
  vector_type x_a;
  vector_type x_b;
  vector_type x;
  x.resize(10,false);
  
  x_a.resize(0,false);
  x_b.resize(10,false);
  
  
  x_b(0) =  1.0;
  x_b(1) =  2.0;
  x_b(2) =  3.0;
  x_b(3) =  4.0;
  x_b(4) =  5.0;
  x_b(5) =  6.0;
  x_b(6) =  7.0;
  x_b(7) =  8.0;
  x_b(8) =  9.0;
  x_b(9) = 10.0;
  
  x.clear();
  
  big::agglomerate_vector( x_a, x_b, new2old, x);
  
  double tol = 0.01;
  
  BOOST_CHECK_CLOSE( double( x( 0 ) ), double( x_b(0) ), tol );
  BOOST_CHECK_CLOSE( double( x( 1 ) ), double( x_b(1) ), tol );
  BOOST_CHECK_CLOSE( double( x( 2 ) ), double( x_b(2) ), tol );
  BOOST_CHECK_CLOSE( double( x( 3 ) ), double( x_b(3) ), tol );
  BOOST_CHECK_CLOSE( double( x( 4 ) ), double( x_b(4) ), tol );
  BOOST_CHECK_CLOSE( double( x( 5 ) ), double( x_b(5) ), tol );
  BOOST_CHECK_CLOSE( double( x( 6 ) ), double( x_b(6) ), tol );
  BOOST_CHECK_CLOSE( double( x( 7 ) ), double( x_b(7) ), tol );
  BOOST_CHECK_CLOSE( double( x( 8 ) ), double( x_b(8) ), tol );
  BOOST_CHECK_CLOSE( double( x( 9 ) ), double( x_b(9) ), tol );
}

BOOST_AUTO_TEST_CASE(test_case3)
{
  typedef ublas::vector<size_t>       idx_vector_type;
  typedef ublas::vector<double>       vector_type;
  
  idx_vector_type new2old;
  idx_vector_type old2new;
  
  old2new.resize(10,false);
  new2old.resize(10,false);
  
  old2new( 0 ) = 0 ;
  old2new( 1 ) = 1 ;
  old2new( 2 ) = 2 ;
  old2new( 3 ) = 3 ;
  old2new( 4 ) = 4 ;
  old2new( 5 ) = 5 ;
  old2new( 6 ) = 6 ;
  old2new( 7 ) = 7 ;
  old2new( 8 ) = 8 ;
  old2new( 9 ) = 9 ;
  
  new2old( 0 ) = 0 ;
  new2old( 1 ) = 1 ;
  new2old( 2 ) = 2 ;
  new2old( 3 ) = 3 ;
  new2old( 4 ) = 4 ;
  new2old( 5 ) = 5 ;
  new2old( 6 ) = 6 ;
  new2old( 7 ) = 7 ;
  new2old( 8 ) = 8 ;
  new2old( 9 ) = 9 ;
  
  vector_type x_a;
  vector_type x_b;
  vector_type x;
  x.resize(10,false);
  
  x_a.resize(10,false);
  x_b.resize(0,false);
  
  
  x_a(0) =  1.0;
  x_a(1) =  2.0;
  x_a(2) =  3.0;
  x_a(3) =  4.0;
  x_a(4) =  5.0;
  x_a(5) =  6.0;
  x_a(6) =  7.0;
  x_a(7) =  8.0;
  x_a(8) =  9.0;
  x_a(9) = 10.0;
  
  x.clear();
  
  big::agglomerate_vector( x_a, x_b, new2old, x);
  
  double tol = 0.01;
  
  BOOST_CHECK_CLOSE( double( x( 0 ) ), double( x_a(0) ), tol );
  BOOST_CHECK_CLOSE( double( x( 1 ) ), double( x_a(1) ), tol );
  BOOST_CHECK_CLOSE( double( x( 2 ) ), double( x_a(2) ), tol );
  BOOST_CHECK_CLOSE( double( x( 3 ) ), double( x_a(3) ), tol );
  BOOST_CHECK_CLOSE( double( x( 4 ) ), double( x_a(4) ), tol );
  BOOST_CHECK_CLOSE( double( x( 5 ) ), double( x_a(5) ), tol );
  BOOST_CHECK_CLOSE( double( x( 6 ) ), double( x_a(6) ), tol );
  BOOST_CHECK_CLOSE( double( x( 7 ) ), double( x_a(7) ), tol );
  BOOST_CHECK_CLOSE( double( x( 8 ) ), double( x_a(8) ), tol );
  BOOST_CHECK_CLOSE( double( x( 9 ) ), double( x_a(9) ), tol );
}


BOOST_AUTO_TEST_SUITE_END();
