#include <big_random.h>
#include <big_return_codes.h>
#include <big_types.h>
#include <min_map/big_compute_jacobian.h>

#define BOOST_AUTO_TEST_MAIN
#include <boost/test/auto_unit_test.hpp>
#include <boost/test/unit_test_suite.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/test_tools.hpp>

BOOST_AUTO_TEST_SUITE(big_compute_full_jacobian);

BOOST_AUTO_TEST_CASE(test_case1)
{
  typedef ublas::vector<size_t>            idx_vector_type;
  typedef ublas::vector<double>            vector_type;
  typedef ublas::compressed_matrix<double> matrix_type;
  typedef big::ValueTraits<double>         value_traits;
  
  typedef double real_type;
  typedef size_t size_type;

  matrix_type A;
  A.resize(10,10,false);

  big::Random<double> value(0.0,1.0);
  for(size_t i=0;i<A.size1();++i)
  {
    for(size_t j=0;j<A.size2();++j)
      A(i,j) = value();
  }

  idx_vector_type bitmask;

  bitmask.resize(10,false);

  bitmask(0) = big::IN_NON_ACTIVE;
  bitmask(1) = big::IN_NON_ACTIVE;
  bitmask(2) = big::IN_ACTIVE;
  bitmask(3) = big::IN_NON_ACTIVE;
  bitmask(4) = big::IN_ACTIVE;
  bitmask(5) = big::IN_NON_ACTIVE;
  bitmask(6) = big::IN_ACTIVE;
  bitmask(7) = big::IN_NON_ACTIVE;
  bitmask(8) = big::IN_NON_ACTIVE;
  bitmask(9) = big::IN_ACTIVE;

  matrix_type J;
  big::compute_jacobian(
      A
    , bitmask        
    , J
    );

  double tol = 0.01;

  // in-active
  BOOST_CHECK_CLOSE( double( J(0,0) ), double( 1.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(0,1) ), double( 0.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(0,2) ), double( 0.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(0,3) ), double( 0.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(0,4) ), double( 0.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(0,5) ), double( 0.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(0,6) ), double( 0.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(0,7) ), double( 0.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(0,8) ), double( 0.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(0,9) ), double( 0.0 ), tol );

  // in-active
  BOOST_CHECK_CLOSE( double( J(1,0) ), double( 0.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(1,1) ), double( 1.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(1,2) ), double( 0.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(1,3) ), double( 0.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(1,4) ), double( 0.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(1,5) ), double( 0.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(1,6) ), double( 0.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(1,7) ), double( 0.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(1,8) ), double( 0.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(1,9) ), double( 0.0 ), tol );

  // active
  BOOST_CHECK_CLOSE( double( J(2,0) ), double( A(2,0) ), tol );
  BOOST_CHECK_CLOSE( double( J(2,1) ), double( A(2,1) ), tol );
  BOOST_CHECK_CLOSE( double( J(2,2) ), double( A(2,2) ), tol );
  BOOST_CHECK_CLOSE( double( J(2,3) ), double( A(2,3) ), tol );
  BOOST_CHECK_CLOSE( double( J(2,4) ), double( A(2,4) ), tol );
  BOOST_CHECK_CLOSE( double( J(2,5) ), double( A(2,5) ), tol );
  BOOST_CHECK_CLOSE( double( J(2,6) ), double( A(2,6) ), tol );
  BOOST_CHECK_CLOSE( double( J(2,7) ), double( A(2,7) ), tol );
  BOOST_CHECK_CLOSE( double( J(2,8) ), double( A(2,8) ), tol );
  BOOST_CHECK_CLOSE( double( J(2,9) ), double( A(2,9) ), tol );

  // in-active
  BOOST_CHECK_CLOSE( double( J(3,0) ), double( 0.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(3,1) ), double( 0.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(3,2) ), double( 0.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(3,3) ), double( 1.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(3,4) ), double( 0.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(3,5) ), double( 0.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(3,6) ), double( 0.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(3,7) ), double( 0.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(3,8) ), double( 0.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(3,9) ), double( 0.0 ), tol );

  // active
  BOOST_CHECK_CLOSE( double( J(4,0) ), double( A(4,0) ), tol );
  BOOST_CHECK_CLOSE( double( J(4,1) ), double( A(4,1) ), tol );
  BOOST_CHECK_CLOSE( double( J(4,2) ), double( A(4,2) ), tol );
  BOOST_CHECK_CLOSE( double( J(4,3) ), double( A(4,3) ), tol );
  BOOST_CHECK_CLOSE( double( J(4,4) ), double( A(4,4) ), tol );
  BOOST_CHECK_CLOSE( double( J(4,5) ), double( A(4,5) ), tol );
  BOOST_CHECK_CLOSE( double( J(4,6) ), double( A(4,6) ), tol );
  BOOST_CHECK_CLOSE( double( J(4,7) ), double( A(4,7) ), tol );
  BOOST_CHECK_CLOSE( double( J(4,8) ), double( A(4,8) ), tol );
  BOOST_CHECK_CLOSE( double( J(4,9) ), double( A(4,9) ), tol );

  // in-active
  BOOST_CHECK_CLOSE( double( J(5,0) ), double( 0.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(5,1) ), double( 0.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(5,2) ), double( 0.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(5,3) ), double( 0.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(5,4) ), double( 0.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(5,5) ), double( 1.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(5,6) ), double( 0.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(5,7) ), double( 0.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(5,8) ), double( 0.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(5,9) ), double( 0.0 ), tol );

  // active
  BOOST_CHECK_CLOSE( double( J(6,0) ), double( A(6,0) ), tol );
  BOOST_CHECK_CLOSE( double( J(6,1) ), double( A(6,1) ), tol );
  BOOST_CHECK_CLOSE( double( J(6,2) ), double( A(6,2) ), tol );
  BOOST_CHECK_CLOSE( double( J(6,3) ), double( A(6,3) ), tol );
  BOOST_CHECK_CLOSE( double( J(6,4) ), double( A(6,4) ), tol );
  BOOST_CHECK_CLOSE( double( J(6,5) ), double( A(6,5) ), tol );
  BOOST_CHECK_CLOSE( double( J(6,6) ), double( A(6,6) ), tol );
  BOOST_CHECK_CLOSE( double( J(6,7) ), double( A(6,7) ), tol );
  BOOST_CHECK_CLOSE( double( J(6,8) ), double( A(6,8) ), tol );
  BOOST_CHECK_CLOSE( double( J(6,9) ), double( A(6,9) ), tol );

  // in-active
  BOOST_CHECK_CLOSE( double( J(7,0) ), double( 0.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(7,1) ), double( 0.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(7,2) ), double( 0.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(7,3) ), double( 0.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(7,4) ), double( 0.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(7,5) ), double( 0.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(7,6) ), double( 0.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(7,7) ), double( 1.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(7,8) ), double( 0.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(7,9) ), double( 0.0 ), tol );

  // in-active
  BOOST_CHECK_CLOSE( double( J(8,0) ), double( 0.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(8,1) ), double( 0.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(8,2) ), double( 0.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(8,3) ), double( 0.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(8,4) ), double( 0.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(8,5) ), double( 0.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(8,6) ), double( 0.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(8,7) ), double( 0.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(8,8) ), double( 1.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(8,9) ), double( 0.0 ), tol );

  // active
  BOOST_CHECK_CLOSE( double( J(9,0) ), double( A(9,0) ), tol );
  BOOST_CHECK_CLOSE( double( J(9,1) ), double( A(9,1) ), tol );
  BOOST_CHECK_CLOSE( double( J(9,2) ), double( A(9,2) ), tol );
  BOOST_CHECK_CLOSE( double( J(9,3) ), double( A(9,3) ), tol );
  BOOST_CHECK_CLOSE( double( J(9,4) ), double( A(9,4) ), tol );
  BOOST_CHECK_CLOSE( double( J(9,5) ), double( A(9,5) ), tol );
  BOOST_CHECK_CLOSE( double( J(9,6) ), double( A(9,6) ), tol );
  BOOST_CHECK_CLOSE( double( J(9,7) ), double( A(9,7) ), tol );
  BOOST_CHECK_CLOSE( double( J(9,8) ), double( A(9,8) ), tol );
  BOOST_CHECK_CLOSE( double( J(9,9) ), double( A(9,9) ), tol );

}

BOOST_AUTO_TEST_CASE(test_case2)
{
  typedef ublas::vector<size_t>            idx_vector_type;
  typedef ublas::vector<double>            vector_type;
  typedef ublas::compressed_matrix<double> matrix_type;
  typedef big::ValueTraits<double>         value_traits;
  
  typedef double real_type;
  typedef size_t size_type;
  
  matrix_type A;
  A.resize(10,10,false);
  
  big::Random<double> value(0.0,1.0);
  for(size_t i=0;i<A.size1();++i)
  {
    for(size_t j=0;j<A.size2();++j)
      A(i,j) = value();
  }
  
  idx_vector_type bitmask;
  
  bitmask.resize(10,false);
  
  bitmask(0) = big::IN_NON_ACTIVE;
  bitmask(1) = big::IN_NON_ACTIVE;
  bitmask(2) = big::IN_NON_ACTIVE;
  bitmask(3) = big::IN_NON_ACTIVE;
  bitmask(4) = big::IN_NON_ACTIVE;
  bitmask(5) = big::IN_NON_ACTIVE;
  bitmask(6) = big::IN_NON_ACTIVE;
  bitmask(7) = big::IN_NON_ACTIVE;
  bitmask(8) = big::IN_NON_ACTIVE;
  bitmask(9) = big::IN_NON_ACTIVE;
  
  matrix_type J;
  big::compute_jacobian(
                        A
                        , bitmask
                        , J
                        );
  
  double tol = 0.01;
  
  // in-active
  BOOST_CHECK_CLOSE( double( J(0,0) ), double( 1.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(0,1) ), double( 0.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(0,2) ), double( 0.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(0,3) ), double( 0.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(0,4) ), double( 0.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(0,5) ), double( 0.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(0,6) ), double( 0.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(0,7) ), double( 0.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(0,8) ), double( 0.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(0,9) ), double( 0.0 ), tol );
  
  // in-active
  BOOST_CHECK_CLOSE( double( J(1,0) ), double( 0.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(1,1) ), double( 1.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(1,2) ), double( 0.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(1,3) ), double( 0.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(1,4) ), double( 0.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(1,5) ), double( 0.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(1,6) ), double( 0.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(1,7) ), double( 0.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(1,8) ), double( 0.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(1,9) ), double( 0.0 ), tol );
  
  // in-active
  BOOST_CHECK_CLOSE( double( J(2,0) ), double( 0.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(2,1) ), double( 0.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(2,2) ), double( 1.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(2,3) ), double( 0.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(2,4) ), double( 0.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(2,5) ), double( 0.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(2,6) ), double( 0.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(2,7) ), double( 0.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(2,8) ), double( 0.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(2,9) ), double( 0.0 ), tol );
  
  // in-active
  BOOST_CHECK_CLOSE( double( J(3,0) ), double( 0.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(3,1) ), double( 0.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(3,2) ), double( 0.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(3,3) ), double( 1.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(3,4) ), double( 0.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(3,5) ), double( 0.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(3,6) ), double( 0.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(3,7) ), double( 0.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(3,8) ), double( 0.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(3,9) ), double( 0.0 ), tol );
  
  // in-active
  BOOST_CHECK_CLOSE( double( J(4,0) ), double( 0.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(4,1) ), double( 0.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(4,2) ), double( 0.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(4,3) ), double( 0.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(4,4) ), double( 1.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(4,5) ), double( 0.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(4,6) ), double( 0.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(4,7) ), double( 0.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(4,8) ), double( 0.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(4,9) ), double( 0.0 ), tol );
  
  // in-active
  BOOST_CHECK_CLOSE( double( J(5,0) ), double( 0.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(5,1) ), double( 0.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(5,2) ), double( 0.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(5,3) ), double( 0.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(5,4) ), double( 0.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(5,5) ), double( 1.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(5,6) ), double( 0.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(5,7) ), double( 0.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(5,8) ), double( 0.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(5,9) ), double( 0.0 ), tol );
  
  // in-active
  BOOST_CHECK_CLOSE( double( J(6,0) ), double( 0.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(6,1) ), double( 0.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(6,2) ), double( 0.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(6,3) ), double( 0.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(6,4) ), double( 0.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(6,5) ), double( 0.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(6,6) ), double( 1.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(6,7) ), double( 0.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(6,8) ), double( 0.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(6,9) ), double( 0.0 ), tol );
  
  // in-active
  BOOST_CHECK_CLOSE( double( J(7,0) ), double( 0.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(7,1) ), double( 0.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(7,2) ), double( 0.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(7,3) ), double( 0.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(7,4) ), double( 0.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(7,5) ), double( 0.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(7,6) ), double( 0.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(7,7) ), double( 1.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(7,8) ), double( 0.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(7,9) ), double( 0.0 ), tol );
  
  // in-active
  BOOST_CHECK_CLOSE( double( J(8,0) ), double( 0.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(8,1) ), double( 0.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(8,2) ), double( 0.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(8,3) ), double( 0.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(8,4) ), double( 0.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(8,5) ), double( 0.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(8,6) ), double( 0.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(8,7) ), double( 0.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(8,8) ), double( 1.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(8,9) ), double( 0.0 ), tol );
  
  // in-active
  BOOST_CHECK_CLOSE( double( J(9,0) ), double( 0.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(9,1) ), double( 0.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(9,2) ), double( 0.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(9,3) ), double( 0.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(9,4) ), double( 0.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(9,5) ), double( 0.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(9,6) ), double( 0.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(9,7) ), double( 0.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(9,8) ), double( 0.0 ), tol );
  BOOST_CHECK_CLOSE( double( J(9,9) ), double( 1.0 ), tol );
  
}

BOOST_AUTO_TEST_CASE(test_case3)
{
  typedef ublas::vector<size_t>            idx_vector_type;
  typedef ublas::vector<double>            vector_type;
  typedef ublas::compressed_matrix<double> matrix_type;
  typedef big::ValueTraits<double>         value_traits;
  
  typedef double real_type;
  typedef size_t size_type;
  
  matrix_type A;
  A.resize(10,10,false);
  
  big::Random<double> value(0.0,1.0);
  for(size_t i=0;i<A.size1();++i)
  {
    for(size_t j=0;j<A.size2();++j)
      A(i,j) = value();
  }
  
  idx_vector_type bitmask;
  
  bitmask.resize(10,false);
  
  bitmask(0) = big::IN_ACTIVE;
  bitmask(1) = big::IN_ACTIVE;
  bitmask(2) = big::IN_ACTIVE;
  bitmask(3) = big::IN_ACTIVE;
  bitmask(4) = big::IN_ACTIVE;
  bitmask(5) = big::IN_ACTIVE;
  bitmask(6) = big::IN_ACTIVE;
  bitmask(7) = big::IN_ACTIVE;
  bitmask(8) = big::IN_ACTIVE;
  bitmask(9) = big::IN_ACTIVE;
  
  matrix_type J;
  big::compute_jacobian(
                        A
                        , bitmask
                        , J
                        );
  
  double tol = 0.01;
  
  // active
  BOOST_CHECK_CLOSE( double( J(0,0) ), double( A(0,0) ), tol );
  BOOST_CHECK_CLOSE( double( J(0,1) ), double( A(0,1) ), tol );
  BOOST_CHECK_CLOSE( double( J(0,2) ), double( A(0,2) ), tol );
  BOOST_CHECK_CLOSE( double( J(0,3) ), double( A(0,3) ), tol );
  BOOST_CHECK_CLOSE( double( J(0,4) ), double( A(0,4) ), tol );
  BOOST_CHECK_CLOSE( double( J(0,5) ), double( A(0,5) ), tol );
  BOOST_CHECK_CLOSE( double( J(0,6) ), double( A(0,6) ), tol );
  BOOST_CHECK_CLOSE( double( J(0,7) ), double( A(0,7) ), tol );
  BOOST_CHECK_CLOSE( double( J(0,8) ), double( A(0,8) ), tol );
  BOOST_CHECK_CLOSE( double( J(0,9) ), double( A(0,9) ), tol );
  
  // active
  BOOST_CHECK_CLOSE( double( J(1,0) ), double( A(1,0) ), tol );
  BOOST_CHECK_CLOSE( double( J(1,1) ), double( A(1,1) ), tol );
  BOOST_CHECK_CLOSE( double( J(1,2) ), double( A(1,2) ), tol );
  BOOST_CHECK_CLOSE( double( J(1,3) ), double( A(1,3) ), tol );
  BOOST_CHECK_CLOSE( double( J(1,4) ), double( A(1,4) ), tol );
  BOOST_CHECK_CLOSE( double( J(1,5) ), double( A(1,5) ), tol );
  BOOST_CHECK_CLOSE( double( J(1,6) ), double( A(1,6) ), tol );
  BOOST_CHECK_CLOSE( double( J(1,7) ), double( A(1,7) ), tol );
  BOOST_CHECK_CLOSE( double( J(1,8) ), double( A(1,8) ), tol );
  BOOST_CHECK_CLOSE( double( J(1,9) ), double( A(1,9) ), tol );
  
  // active
  BOOST_CHECK_CLOSE( double( J(2,0) ), double( A(2,0) ), tol );
  BOOST_CHECK_CLOSE( double( J(2,1) ), double( A(2,1) ), tol );
  BOOST_CHECK_CLOSE( double( J(2,2) ), double( A(2,2) ), tol );
  BOOST_CHECK_CLOSE( double( J(2,3) ), double( A(2,3) ), tol );
  BOOST_CHECK_CLOSE( double( J(2,4) ), double( A(2,4) ), tol );
  BOOST_CHECK_CLOSE( double( J(2,5) ), double( A(2,5) ), tol );
  BOOST_CHECK_CLOSE( double( J(2,6) ), double( A(2,6) ), tol );
  BOOST_CHECK_CLOSE( double( J(2,7) ), double( A(2,7) ), tol );
  BOOST_CHECK_CLOSE( double( J(2,8) ), double( A(2,8) ), tol );
  BOOST_CHECK_CLOSE( double( J(2,9) ), double( A(2,9) ), tol );
  
  // active
  BOOST_CHECK_CLOSE( double( J(3,0) ), double( A(3,0) ), tol );
  BOOST_CHECK_CLOSE( double( J(3,1) ), double( A(3,1) ), tol );
  BOOST_CHECK_CLOSE( double( J(3,2) ), double( A(3,2) ), tol );
  BOOST_CHECK_CLOSE( double( J(3,3) ), double( A(3,3) ), tol );
  BOOST_CHECK_CLOSE( double( J(3,4) ), double( A(3,4) ), tol );
  BOOST_CHECK_CLOSE( double( J(3,5) ), double( A(3,5) ), tol );
  BOOST_CHECK_CLOSE( double( J(3,6) ), double( A(3,6) ), tol );
  BOOST_CHECK_CLOSE( double( J(3,7) ), double( A(3,7) ), tol );
  BOOST_CHECK_CLOSE( double( J(3,8) ), double( A(3,8) ), tol );
  BOOST_CHECK_CLOSE( double( J(3,9) ), double( A(3,9) ), tol );
  
  // active
  BOOST_CHECK_CLOSE( double( J(4,0) ), double( A(4,0) ), tol );
  BOOST_CHECK_CLOSE( double( J(4,1) ), double( A(4,1) ), tol );
  BOOST_CHECK_CLOSE( double( J(4,2) ), double( A(4,2) ), tol );
  BOOST_CHECK_CLOSE( double( J(4,3) ), double( A(4,3) ), tol );
  BOOST_CHECK_CLOSE( double( J(4,4) ), double( A(4,4) ), tol );
  BOOST_CHECK_CLOSE( double( J(4,5) ), double( A(4,5) ), tol );
  BOOST_CHECK_CLOSE( double( J(4,6) ), double( A(4,6) ), tol );
  BOOST_CHECK_CLOSE( double( J(4,7) ), double( A(4,7) ), tol );
  BOOST_CHECK_CLOSE( double( J(4,8) ), double( A(4,8) ), tol );
  BOOST_CHECK_CLOSE( double( J(4,9) ), double( A(4,9) ), tol );
  
  // active
  BOOST_CHECK_CLOSE( double( J(5,0) ), double( A(5,0) ), tol );
  BOOST_CHECK_CLOSE( double( J(5,1) ), double( A(5,1) ), tol );
  BOOST_CHECK_CLOSE( double( J(5,2) ), double( A(5,2) ), tol );
  BOOST_CHECK_CLOSE( double( J(5,3) ), double( A(5,3) ), tol );
  BOOST_CHECK_CLOSE( double( J(5,4) ), double( A(5,4) ), tol );
  BOOST_CHECK_CLOSE( double( J(5,5) ), double( A(5,5) ), tol );
  BOOST_CHECK_CLOSE( double( J(5,6) ), double( A(5,6) ), tol );
  BOOST_CHECK_CLOSE( double( J(5,7) ), double( A(5,7) ), tol );
  BOOST_CHECK_CLOSE( double( J(5,8) ), double( A(5,8) ), tol );
  BOOST_CHECK_CLOSE( double( J(5,9) ), double( A(5,9) ), tol );
  
  // active
  BOOST_CHECK_CLOSE( double( J(6,0) ), double( A(6,0) ), tol );
  BOOST_CHECK_CLOSE( double( J(6,1) ), double( A(6,1) ), tol );
  BOOST_CHECK_CLOSE( double( J(6,2) ), double( A(6,2) ), tol );
  BOOST_CHECK_CLOSE( double( J(6,3) ), double( A(6,3) ), tol );
  BOOST_CHECK_CLOSE( double( J(6,4) ), double( A(6,4) ), tol );
  BOOST_CHECK_CLOSE( double( J(6,5) ), double( A(6,5) ), tol );
  BOOST_CHECK_CLOSE( double( J(6,6) ), double( A(6,6) ), tol );
  BOOST_CHECK_CLOSE( double( J(6,7) ), double( A(6,7) ), tol );
  BOOST_CHECK_CLOSE( double( J(6,8) ), double( A(6,8) ), tol );
  BOOST_CHECK_CLOSE( double( J(6,9) ), double( A(6,9) ), tol );
  
  // active
  BOOST_CHECK_CLOSE( double( J(7,0) ), double( A(7,0) ), tol );
  BOOST_CHECK_CLOSE( double( J(7,1) ), double( A(7,1) ), tol );
  BOOST_CHECK_CLOSE( double( J(7,2) ), double( A(7,2) ), tol );
  BOOST_CHECK_CLOSE( double( J(7,3) ), double( A(7,3) ), tol );
  BOOST_CHECK_CLOSE( double( J(7,4) ), double( A(7,4) ), tol );
  BOOST_CHECK_CLOSE( double( J(7,5) ), double( A(7,5) ), tol );
  BOOST_CHECK_CLOSE( double( J(7,6) ), double( A(7,6) ), tol );
  BOOST_CHECK_CLOSE( double( J(7,7) ), double( A(7,7) ), tol );
  BOOST_CHECK_CLOSE( double( J(7,8) ), double( A(7,8) ), tol );
  BOOST_CHECK_CLOSE( double( J(7,9) ), double( A(7,9) ), tol );
  
  // active
  BOOST_CHECK_CLOSE( double( J(8,0) ), double( A(8,0) ), tol );
  BOOST_CHECK_CLOSE( double( J(8,1) ), double( A(8,1) ), tol );
  BOOST_CHECK_CLOSE( double( J(8,2) ), double( A(8,2) ), tol );
  BOOST_CHECK_CLOSE( double( J(8,3) ), double( A(8,3) ), tol );
  BOOST_CHECK_CLOSE( double( J(8,4) ), double( A(8,4) ), tol );
  BOOST_CHECK_CLOSE( double( J(8,5) ), double( A(8,5) ), tol );
  BOOST_CHECK_CLOSE( double( J(8,6) ), double( A(8,6) ), tol );
  BOOST_CHECK_CLOSE( double( J(8,7) ), double( A(8,7) ), tol );
  BOOST_CHECK_CLOSE( double( J(8,8) ), double( A(8,8) ), tol );
  BOOST_CHECK_CLOSE( double( J(8,9) ), double( A(8,9) ), tol );
  
  // active
  BOOST_CHECK_CLOSE( double( J(9,0) ), double( A(9,0) ), tol );
  BOOST_CHECK_CLOSE( double( J(9,1) ), double( A(9,1) ), tol );
  BOOST_CHECK_CLOSE( double( J(9,2) ), double( A(9,2) ), tol );
  BOOST_CHECK_CLOSE( double( J(9,3) ), double( A(9,3) ), tol );
  BOOST_CHECK_CLOSE( double( J(9,4) ), double( A(9,4) ), tol );
  BOOST_CHECK_CLOSE( double( J(9,5) ), double( A(9,5) ), tol );
  BOOST_CHECK_CLOSE( double( J(9,6) ), double( A(9,6) ), tol );
  BOOST_CHECK_CLOSE( double( J(9,7) ), double( A(9,7) ), tol );
  BOOST_CHECK_CLOSE( double( J(9,8) ), double( A(9,8) ), tol );
  BOOST_CHECK_CLOSE( double( J(9,9) ), double( A(9,9) ), tol );
}

BOOST_AUTO_TEST_SUITE_END();
