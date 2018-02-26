#include <big_types.h>
#include <big_return_codes.h>
#include <min_map/big_compute_index_reordering.h>

#define BOOST_AUTO_TEST_MAIN
#include <boost/test/auto_unit_test.hpp>
#include <boost/test/unit_test_suite.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/test_tools.hpp>

BOOST_AUTO_TEST_SUITE(big_compute_index_reordering);

BOOST_AUTO_TEST_CASE(test_case1)
{
  typedef ublas::vector<size_t>   idx_vector_type;

  idx_vector_type bitmask;

  bitmask.resize(10,false);

  bitmask(2) = big::IN_ACTIVE;            // 0
  bitmask(4) = big::IN_ACTIVE;            // 1
  bitmask(6) = big::IN_ACTIVE;            // 2
  bitmask(9) = big::IN_ACTIVE;            // 3
  bitmask(0) = big::IN_NON_ACTIVE;        // 4
  bitmask(1) = big::IN_NON_ACTIVE;        // 5
  bitmask(3) = big::IN_NON_ACTIVE;        // 6
  bitmask(5) = big::IN_NON_ACTIVE;        // 7
  bitmask(7) = big::IN_NON_ACTIVE;        // 8
  bitmask(8) = big::IN_NON_ACTIVE;        // 9
  
  idx_vector_type old2new;
  idx_vector_type new2old;

  big::compute_index_reordering( bitmask, old2new, new2old );

  BOOST_CHECK( old2new( 2 ) == 0 );
  BOOST_CHECK( old2new( 4 ) == 1 );
  BOOST_CHECK( old2new( 6 ) == 2 );
  BOOST_CHECK( old2new( 9 ) == 3 );
  BOOST_CHECK( old2new( 0 ) == 4 );
  BOOST_CHECK( old2new( 1 ) == 5 );
  BOOST_CHECK( old2new( 3 ) == 6 );
  BOOST_CHECK( old2new( 5 ) == 7 );
  BOOST_CHECK( old2new( 7 ) == 8 );
  BOOST_CHECK( old2new( 8 ) == 9 );
  
  BOOST_CHECK( new2old( 0 ) == 2 );
  BOOST_CHECK( new2old( 1 ) == 4 );
  BOOST_CHECK( new2old( 2 ) == 6 );
  BOOST_CHECK( new2old( 3 ) == 9 );
  BOOST_CHECK( new2old( 4 ) == 0 );
  BOOST_CHECK( new2old( 5 ) == 1 );
  BOOST_CHECK( new2old( 6 ) == 3 );
  BOOST_CHECK( new2old( 7 ) == 5 );
  BOOST_CHECK( new2old( 8 ) == 7 );
  BOOST_CHECK( new2old( 9 ) == 8 );

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
  
  bitmask(1) = big::IN_ACTIVE;
  bitmask(3) = big::IN_ACTIVE;
  bitmask(8) = big::IN_ACTIVE;
  bitmask(0) = big::IN_ACTIVE;
  bitmask(5) = big::IN_ACTIVE;
  bitmask(7) = big::IN_ACTIVE;
  
  idx_vector_type old2new;
  idx_vector_type new2old;
  
  big::compute_index_reordering( bitmask, old2new, new2old );
  
  BOOST_CHECK( old2new( 0 ) == 0 );
  BOOST_CHECK( old2new( 1 ) == 1 );
  BOOST_CHECK( old2new( 2 ) == 2 );
  BOOST_CHECK( old2new( 3 ) == 3 );
  BOOST_CHECK( old2new( 4 ) == 4 );
  BOOST_CHECK( old2new( 5 ) == 5 );
  BOOST_CHECK( old2new( 6 ) == 6 );
  BOOST_CHECK( old2new( 7 ) == 7 );
  BOOST_CHECK( old2new( 8 ) == 8 );
  BOOST_CHECK( old2new( 9 ) == 9 );
  
  
  BOOST_CHECK( new2old( 0 ) == 0 );
  BOOST_CHECK( new2old( 1 ) == 1 );
  BOOST_CHECK( new2old( 2 ) == 2 );
  BOOST_CHECK( new2old( 3 ) == 3 );
  BOOST_CHECK( new2old( 4 ) == 4 );
  BOOST_CHECK( new2old( 5 ) == 5 );
  BOOST_CHECK( new2old( 6 ) == 6 );
  BOOST_CHECK( new2old( 7 ) == 7 );
  BOOST_CHECK( new2old( 8 ) == 8 );
  BOOST_CHECK( new2old( 9 ) == 9 );
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
  
  bitmask(1) = big::IN_NON_ACTIVE;
  bitmask(3) = big::IN_NON_ACTIVE;
  bitmask(8) = big::IN_NON_ACTIVE;
  bitmask(0) = big::IN_NON_ACTIVE;
  bitmask(5) = big::IN_NON_ACTIVE;
  bitmask(7) = big::IN_NON_ACTIVE;
  
  idx_vector_type old2new;
  idx_vector_type new2old;
  
  big::compute_index_reordering( bitmask, old2new, new2old );
  
  BOOST_CHECK( old2new( 0 ) == 0 );
  BOOST_CHECK( old2new( 1 ) == 1 );
  BOOST_CHECK( old2new( 2 ) == 2 );
  BOOST_CHECK( old2new( 3 ) == 3 );
  BOOST_CHECK( old2new( 4 ) == 4 );
  BOOST_CHECK( old2new( 5 ) == 5 );
  BOOST_CHECK( old2new( 6 ) == 6 );
  BOOST_CHECK( old2new( 7 ) == 7 );
  BOOST_CHECK( old2new( 8 ) == 8 );
  BOOST_CHECK( old2new( 9 ) == 9 );
  
  
  BOOST_CHECK( new2old( 0 ) == 0 );
  BOOST_CHECK( new2old( 1 ) == 1 );
  BOOST_CHECK( new2old( 2 ) == 2 );
  BOOST_CHECK( new2old( 3 ) == 3 );
  BOOST_CHECK( new2old( 4 ) == 4 );
  BOOST_CHECK( new2old( 5 ) == 5 );
  BOOST_CHECK( new2old( 6 ) == 6 );
  BOOST_CHECK( new2old( 7 ) == 7 );
  BOOST_CHECK( new2old( 8 ) == 8 );
  BOOST_CHECK( new2old( 9 ) == 9 );
}

BOOST_AUTO_TEST_SUITE_END();
