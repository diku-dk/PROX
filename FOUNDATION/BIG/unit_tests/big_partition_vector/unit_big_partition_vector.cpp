#include <big_types.h>
#include <big_return_codes.h>

#include <min_map/big_compute_index_reordering.h>
#include <min_map/big_partition_vector.h>

#define BOOST_AUTO_TEST_MAIN
#include <boost/test/auto_unit_test.hpp>
#include <boost/test/unit_test_suite.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/test_tools.hpp>

BOOST_AUTO_TEST_SUITE(big_partition_vector);

BOOST_AUTO_TEST_CASE(test_case1)
{
  typedef ublas::vector<size_t>   idx_vector_type;

  idx_vector_type bitmask;

  bitmask.resize(10,false);

  bitmask(2) = big::IN_ACTIVE;
  bitmask(4) = big::IN_ACTIVE;
  bitmask(6) = big::IN_ACTIVE;
  bitmask(9) = big::IN_ACTIVE;
  
  bitmask(0) = big::IN_NON_ACTIVE;
  bitmask(1) = big::IN_NON_ACTIVE;
  bitmask(3) = big::IN_NON_ACTIVE;
  bitmask(5) = big::IN_NON_ACTIVE;
  bitmask(7) = big::IN_NON_ACTIVE;
  bitmask(8) = big::IN_NON_ACTIVE;

  idx_vector_type old2new;
  idx_vector_type new2old;

  big::compute_index_reordering( bitmask, old2new, new2old );

  ublas::vector<double> rhs;
  rhs.resize(10,false);

  rhs(0) = 1.0;
  rhs(1) = 2.0;
  rhs(2) = 3.0;
  rhs(3) = 4.0;
  rhs(4) = 5.0;
  rhs(5) = 6.0;
  rhs(6) = 7.0;
  rhs(7) = 8.0;
  rhs(8) = 9.0;
  rhs(9) = 10.0;

  ublas::vector<double> rhs_a;
  ublas::vector<double> rhs_b;
  big::partition_vector( rhs, bitmask, old2new, 4, 6, rhs_a, rhs_b );

  BOOST_CHECK( rhs_a.size() == 4);
  BOOST_CHECK( rhs_b.size() == 6);

  double tol =0.01;
  BOOST_CHECK_CLOSE( double( rhs_a( 0 ) ), double( rhs( new2old(0) ) ), tol );
  BOOST_CHECK_CLOSE( double( rhs_a( 1 ) ), double( rhs( new2old(1) ) ), tol );
  BOOST_CHECK_CLOSE( double( rhs_a( 2 ) ), double( rhs( new2old(2) ) ), tol );
  BOOST_CHECK_CLOSE( double( rhs_a( 3 ) ), double( rhs( new2old(3) ) ), tol );
  BOOST_CHECK_CLOSE( double( rhs_b( 0 ) ), double( rhs( new2old(0+4) ) ), tol );
  BOOST_CHECK_CLOSE( double( rhs_b( 1 ) ), double( rhs( new2old(1+4) ) ), tol );
  BOOST_CHECK_CLOSE( double( rhs_b( 2 ) ), double( rhs( new2old(2+4) ) ), tol );
  BOOST_CHECK_CLOSE( double( rhs_b( 3 ) ), double( rhs( new2old(3+4) ) ), tol );
  BOOST_CHECK_CLOSE( double( rhs_b( 4 ) ), double( rhs( new2old(4+4) ) ), tol );
  BOOST_CHECK_CLOSE( double( rhs_b( 5 ) ), double( rhs( new2old(5+4) ) ), tol );
}


BOOST_AUTO_TEST_CASE(test_case2)
{
  typedef ublas::vector<size_t>   idx_vector_type;
  
  idx_vector_type bitmask;
  
  bitmask.resize(10,false);
  
  bitmask(2) = big::IN_ACTIVE;
  bitmask(4) = big::IN_ACTIVE;
  bitmask(6) = big::IN_ACTIVE;
  bitmask(9) = big::IN_ACTIVE;
  bitmask(0) = big::IN_ACTIVE;
  bitmask(1) = big::IN_ACTIVE;
  bitmask(3) = big::IN_ACTIVE;
  bitmask(5) = big::IN_ACTIVE;
  bitmask(7) = big::IN_ACTIVE;
  bitmask(8) = big::IN_ACTIVE;
  
  idx_vector_type old2new;
  idx_vector_type new2old;
  
  big::compute_index_reordering( bitmask, old2new, new2old );
  
  ublas::vector<double> rhs;
  rhs.resize(10,false);
  
  rhs(0) = 1.0;
  rhs(1) = 2.0;
  rhs(2) = 3.0;
  rhs(3) = 4.0;
  rhs(4) = 5.0;
  rhs(5) = 6.0;
  rhs(6) = 7.0;
  rhs(7) = 8.0;
  rhs(8) = 9.0;
  rhs(9) = 10.0;
  
  ublas::vector<double> rhs_a;
  ublas::vector<double> rhs_b;
  big::partition_vector( rhs, bitmask, old2new, 10, 0, rhs_a, rhs_b );
  
  BOOST_CHECK( rhs_a.size() == 10);
  BOOST_CHECK( rhs_b.size() == 0);
  
  double tol =0.01;
  BOOST_CHECK_CLOSE( double( rhs_a( 0 ) ), double( rhs( new2old(0) ) ), tol );
  BOOST_CHECK_CLOSE( double( rhs_a( 1 ) ), double( rhs( new2old(1) ) ), tol );
  BOOST_CHECK_CLOSE( double( rhs_a( 2 ) ), double( rhs( new2old(2) ) ), tol );
  BOOST_CHECK_CLOSE( double( rhs_a( 3 ) ), double( rhs( new2old(3) ) ), tol );
  BOOST_CHECK_CLOSE( double( rhs_a( 4 ) ), double( rhs( new2old(4) ) ), tol );
  BOOST_CHECK_CLOSE( double( rhs_a( 5 ) ), double( rhs( new2old(5) ) ), tol );
  BOOST_CHECK_CLOSE( double( rhs_a( 6 ) ), double( rhs( new2old(6) ) ), tol );
  BOOST_CHECK_CLOSE( double( rhs_a( 7 ) ), double( rhs( new2old(7) ) ), tol );
  BOOST_CHECK_CLOSE( double( rhs_a( 8 ) ), double( rhs( new2old(8) ) ), tol );
  BOOST_CHECK_CLOSE( double( rhs_a( 9 ) ), double( rhs( new2old(9) ) ), tol );  
}


BOOST_AUTO_TEST_CASE(test_case3)
{
  typedef ublas::vector<size_t>   idx_vector_type;
  
  idx_vector_type bitmask;
  
  bitmask.resize(10,false);
  
  bitmask(2) = big::IN_NON_ACTIVE;
  bitmask(4) = big::IN_NON_ACTIVE;
  bitmask(6) = big::IN_NON_ACTIVE;
  bitmask(9) = big::IN_NON_ACTIVE;
  bitmask(0) = big::IN_NON_ACTIVE;
  bitmask(1) = big::IN_NON_ACTIVE;
  bitmask(3) = big::IN_NON_ACTIVE;
  bitmask(5) = big::IN_NON_ACTIVE;
  bitmask(7) = big::IN_NON_ACTIVE;
  bitmask(8) = big::IN_NON_ACTIVE;
  
  idx_vector_type old2new;
  idx_vector_type new2old;
  
  big::compute_index_reordering( bitmask, old2new, new2old );
  
  ublas::vector<double> rhs;
  rhs.resize(10,false);
  
  rhs(0) = 1.0;
  rhs(1) = 2.0;
  rhs(2) = 3.0;
  rhs(3) = 4.0;
  rhs(4) = 5.0;
  rhs(5) = 6.0;
  rhs(6) = 7.0;
  rhs(7) = 8.0;
  rhs(8) = 9.0;
  rhs(9) = 10.0;
  
  ublas::vector<double> rhs_a;
  ublas::vector<double> rhs_b;
  big::partition_vector( rhs, bitmask, old2new, 0, 10, rhs_a, rhs_b );
  
  BOOST_CHECK( rhs_a.size() == 0);
  BOOST_CHECK( rhs_b.size() == 10);
  
  double tol =0.01;
  BOOST_CHECK_CLOSE( double( rhs_b( 0 ) ), double( rhs( new2old(0) ) ), tol );
  BOOST_CHECK_CLOSE( double( rhs_b( 1 ) ), double( rhs( new2old(1) ) ), tol );
  BOOST_CHECK_CLOSE( double( rhs_b( 2 ) ), double( rhs( new2old(2) ) ), tol );
  BOOST_CHECK_CLOSE( double( rhs_b( 3 ) ), double( rhs( new2old(3) ) ), tol );
  BOOST_CHECK_CLOSE( double( rhs_b( 4 ) ), double( rhs( new2old(4) ) ), tol );
  BOOST_CHECK_CLOSE( double( rhs_b( 5 ) ), double( rhs( new2old(5) ) ), tol );
  BOOST_CHECK_CLOSE( double( rhs_b( 6 ) ), double( rhs( new2old(6) ) ), tol );
  BOOST_CHECK_CLOSE( double( rhs_b( 7 ) ), double( rhs( new2old(7) ) ), tol );
  BOOST_CHECK_CLOSE( double( rhs_b( 8 ) ), double( rhs( new2old(8) ) ), tol );
  BOOST_CHECK_CLOSE( double( rhs_b( 9 ) ), double( rhs( new2old(9) ) ), tol );
}

BOOST_AUTO_TEST_SUITE_END();
