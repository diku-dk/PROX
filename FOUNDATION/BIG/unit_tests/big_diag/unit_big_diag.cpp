#include <big_types.h>
#include <big_generate_random.h>
#include <big_diag.h>

#define BOOST_AUTO_TEST_MAIN
#include <boost/test/auto_unit_test.hpp>
#include <boost/test/unit_test_suite.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/test_tools.hpp>

BOOST_AUTO_TEST_SUITE(big_diag);

BOOST_AUTO_TEST_CASE(random_test_case)
{
  typedef ublas::compressed_matrix<double> matrix_type;
  typedef ublas::vector<double>            vector_type;
  
  
  for(size_t tst=0u;tst<5u;++tst)
  {
    matrix_type D;
    vector_type v;
    
    big::generate_random(10, v);
    big::diag(v,D);
    
    BOOST_CHECK( D.size1() == 10 );
    BOOST_CHECK( D.size2() == 10 );
    
    for(size_t i=0u;i<10u;++i)
    {
      for(size_t j=0u;j<10u;++j)
      {
        if(i==j)
          BOOST_CHECK( D(i,i) == v(i) );
        else
          BOOST_CHECK( D(i,j) == 0.0 );
      }
    }
  }
}

BOOST_AUTO_TEST_SUITE_END();
