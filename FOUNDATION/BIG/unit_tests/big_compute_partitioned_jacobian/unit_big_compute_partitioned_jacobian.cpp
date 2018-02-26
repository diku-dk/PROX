#include <big_random.h>
#include <big_types.h>
#include <big_return_codes.h>
#include <min_map/big_compute_index_reordering.h>
#include <min_map/big_compute_partitioned_jacobian.h>

#define BOOST_AUTO_TEST_MAIN
#include <boost/test/auto_unit_test.hpp>
#include <boost/test/unit_test_suite.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/test_tools.hpp>

BOOST_AUTO_TEST_SUITE(big_compute_partitioned_jacobian);

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

  idx_vector_type old2new;
  idx_vector_type new2old;

  big::compute_index_reordering( bitmask, old2new, new2old );

  matrix_type A_aa;
  matrix_type A_ab;

  big::compute_partitioned_jacobian(
    A
    ,  bitmask
    , old2new
    ,  4
    ,  6
    ,  A_aa
    ,  A_ab
    );

  double tol = 0.01;

  BOOST_CHECK( A_aa.size1()==4 );
  BOOST_CHECK( A_aa.size2()==4 );
  for(size_type i = 0;i<4;++i)
  {
    for(size_type j = 0;j<4;++j)
    {
      BOOST_CHECK_CLOSE( double( A_aa(i,j) ), double( A(new2old(i),new2old(j)) ), tol );
    }
  }

  BOOST_CHECK( A_ab.size1()==4 );
  BOOST_CHECK( A_ab.size2()==6 );
  for(size_type i = 0;i<4;++i)
  {
    for(size_type j = 0;j<6;++j)
    {
      BOOST_CHECK_CLOSE( double( A_ab(i,j) ), double( A(new2old(i),new2old(j+4)) ), tol );
    }
  }

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
  
  idx_vector_type old2new;
  idx_vector_type new2old;
  
  big::compute_index_reordering( bitmask, old2new, new2old );
  
  matrix_type A_aa;
  matrix_type A_ab;
  
  big::compute_partitioned_jacobian(
                                    A
                                    , bitmask
                                    , old2new
                                    , 10
                                    , 0
                                    , A_aa
                                    , A_ab
                                    );
  
  double tol = 0.01;
  
  BOOST_CHECK( A_aa.size1()==10 );
  BOOST_CHECK( A_aa.size2()==10 );
  for(size_type i = 0;i<10;++i)
  {
    for(size_type j = 0;j<10;++j)
    {
      BOOST_CHECK_CLOSE( double( A_aa(i,j) ), double( A(new2old(i),new2old(j)) ), tol );
    }
  }
  
  BOOST_CHECK( A_ab.size1()==10 );
  BOOST_CHECK( A_ab.size2()==0 );
  
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
  
  idx_vector_type old2new;
  idx_vector_type new2old;
  
  big::compute_index_reordering( bitmask, old2new, new2old );
  
  matrix_type A_aa;
  matrix_type A_ab;
  
  big::compute_partitioned_jacobian(
                                    A
                                    , bitmask
                                    , old2new
                                    , 0
                                    , 10
                                    , A_aa
                                    , A_ab
                                    );
    
  BOOST_CHECK( A_aa.size1()==0 );
  BOOST_CHECK( A_aa.size2()==0 );
  
  BOOST_CHECK( A_ab.size1()==0 );
  BOOST_CHECK( A_ab.size2()==10 );
  
}


BOOST_AUTO_TEST_SUITE_END();
