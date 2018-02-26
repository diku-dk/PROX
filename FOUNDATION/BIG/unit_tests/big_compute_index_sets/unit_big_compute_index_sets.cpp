#include <big_types.h>
#include <big_return_codes.h>
#include <min_map/big_compute_index_sets.h>

#define BOOST_AUTO_TEST_MAIN
#include <boost/test/auto_unit_test.hpp>
#include <boost/test/unit_test_suite.hpp>
#include <boost/test/test_tools.hpp>

BOOST_AUTO_TEST_SUITE(big_compute_index_sets);

BOOST_AUTO_TEST_CASE(test_case1)
{
  typedef double                    real_type;
  typedef ublas::vector<real_type>  vector_type;
  typedef ublas::vector<size_t>     idx_vector_type;

  vector_type x;
  vector_type y;

  x.resize(10,false);
  y.resize(10,false);

  x(0) =  0.0;   y(0) =  0.0;  // inactive
  x(1) = -1.0;   y(1) =  1.0;  // inactive
  x(2) =  1.0;   y(2) = -1.0;  // active
  x(3) = -2.0;   y(3) = -1.0;  // inactive
  x(4) =  2.0;   y(4) =  1.0;  // active
  x(5) = -3.0;   y(5) =  0.0;  // inactive
  x(6) =  3.0;   y(6) =  0.0;  // active
  x(7) = -0.5;   y(7) =  0.5;  // inactive
  x(8) =  0.5;   y(8) = -0.5;  // active
  x(9) =  0.25;  y(9) = 0.25;  // inactive

  idx_vector_type bitmask;
  size_t cnt_active;
  size_t cnt_inactive;
  big::compute_index_sets( y, x, bitmask, cnt_active, cnt_inactive );

  BOOST_CHECK( cnt_active == 4 );
  BOOST_CHECK( cnt_inactive == 6 );

  BOOST_CHECK( bitmask(0) == big::IN_NON_ACTIVE );
  BOOST_CHECK( bitmask(1) == big::IN_NON_ACTIVE );
  BOOST_CHECK( bitmask(2) == big::IN_ACTIVE );
  BOOST_CHECK( bitmask(3) == big::IN_NON_ACTIVE );
  BOOST_CHECK( bitmask(4) == big::IN_ACTIVE );
  BOOST_CHECK( bitmask(5) == big::IN_NON_ACTIVE );
  BOOST_CHECK( bitmask(6) == big::IN_ACTIVE );
  BOOST_CHECK( bitmask(7) == big::IN_NON_ACTIVE );
  BOOST_CHECK( bitmask(8) == big::IN_ACTIVE );
  BOOST_CHECK( bitmask(9) == big::IN_NON_ACTIVE );
}


BOOST_AUTO_TEST_CASE(test_case2)
{
  typedef double                    real_type;
  typedef ublas::vector<real_type>  vector_type;
  typedef ublas::vector<size_t>     idx_vector_type;
  
  vector_type x;
  vector_type y;
  
  x.resize(10,false);
  y.resize(10,false);
  
  x(0) =  0.0;   y(0) =  0.0;  // inactive
  x(1) = -1.0;   y(1) =  1.0;  // inactive
  x(2) = -2.0;   y(2) =  2.0;  // inactive
  x(3) = -3.0;   y(3) =  3.0;  // inactive
  x(4) = -4.0;   y(4) =  4.0;  // inactive
  x(5) = -5.0;   y(5) =  5.0;  // inactive
  x(6) = -6.0;   y(6) =  6.0;  // inactive
  x(7) = -7.0;   y(7) =  7.0;  // inactive
  x(8) = -8.0;   y(8) =  8.0;  // inactive
  x(9) = -9.0;   y(9) =  9.0;  // inactive
  
  idx_vector_type bitmask;
  size_t cnt_active;
  size_t cnt_inactive;
  big::compute_index_sets( y, x, bitmask, cnt_active, cnt_inactive );
  
  BOOST_CHECK( cnt_active == 0 );
  BOOST_CHECK( cnt_inactive == 10 );
  
  BOOST_CHECK( bitmask(0) == big::IN_NON_ACTIVE );
  BOOST_CHECK( bitmask(1) == big::IN_NON_ACTIVE );
  BOOST_CHECK( bitmask(2) == big::IN_NON_ACTIVE );
  BOOST_CHECK( bitmask(3) == big::IN_NON_ACTIVE );
  BOOST_CHECK( bitmask(4) == big::IN_NON_ACTIVE );
  BOOST_CHECK( bitmask(5) == big::IN_NON_ACTIVE );
  BOOST_CHECK( bitmask(6) == big::IN_NON_ACTIVE );
  BOOST_CHECK( bitmask(7) == big::IN_NON_ACTIVE );
  BOOST_CHECK( bitmask(8) == big::IN_NON_ACTIVE );
  BOOST_CHECK( bitmask(9) == big::IN_NON_ACTIVE );
}

BOOST_AUTO_TEST_CASE(test_case3)
{
  typedef double                    real_type;
  typedef ublas::vector<real_type>  vector_type;
  typedef ublas::vector<size_t>     idx_vector_type;
  
  vector_type x;
  vector_type y;
  
  x.resize(10,false);
  y.resize(10,false);
  
  x(0) =  0.1;   y(0) =  0.0;  // active
  x(1) =  1.0;   y(1) = -1.0;  // active
  x(2) =  2.0;   y(2) = -2.0;  // active
  x(3) =  3.0;   y(3) = -3.0;  // active
  x(4) =  4.0;   y(4) = -4.0;  // active
  x(5) =  5.0;   y(5) = -5.0;  // active
  x(6) =  6.0;   y(6) = -6.0;  // active
  x(7) =  7.0;   y(7) = -7.0;  // active
  x(8) =  8.0;   y(8) = -8.0;  // active
  x(9) =  9.0;   y(9) = -9.0;  // active
  
  idx_vector_type bitmask;
  size_t cnt_active;
  size_t cnt_inactive;
  big::compute_index_sets( y, x, bitmask, cnt_active, cnt_inactive );
  
  BOOST_CHECK( cnt_active == 10 );
  BOOST_CHECK( cnt_inactive == 0 );
  
  BOOST_CHECK( bitmask(0) == big::IN_ACTIVE );
  BOOST_CHECK( bitmask(1) == big::IN_ACTIVE );
  BOOST_CHECK( bitmask(2) == big::IN_ACTIVE );
  BOOST_CHECK( bitmask(3) == big::IN_ACTIVE );
  BOOST_CHECK( bitmask(4) == big::IN_ACTIVE );
  BOOST_CHECK( bitmask(5) == big::IN_ACTIVE );
  BOOST_CHECK( bitmask(6) == big::IN_ACTIVE );
  BOOST_CHECK( bitmask(7) == big::IN_ACTIVE );
  BOOST_CHECK( bitmask(8) == big::IN_ACTIVE );
  BOOST_CHECK( bitmask(9) == big::IN_ACTIVE );
}

BOOST_AUTO_TEST_SUITE_END();
