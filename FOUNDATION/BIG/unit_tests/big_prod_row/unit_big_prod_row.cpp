#include <big_random.h>
#include <big_prod_row.h>

#define BOOST_AUTO_TEST_MAIN
#include <boost/test/auto_unit_test.hpp>
#include <boost/test/unit_test_suite.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/test_tools.hpp>

BOOST_AUTO_TEST_SUITE(big_prod_row);

BOOST_AUTO_TEST_CASE(random_test_case)
{
  typedef ublas::compressed_matrix<double> matrix_type;
  typedef ublas::vector<double>            vector_type;

  matrix_type A;
  vector_type x;
  vector_type y;
  vector_type tst;

  A.resize(10,10,false);
  x.resize(10,false);
  y.resize(10,false);
  tst.resize(10,false);

  big::Random<double> value(0.0,1.0);

  for(size_t iteration = 0;iteration<100;++iteration)
  {
    for(size_t i=0;i<A.size1();++i)
    {
      x(i) = value();
      y(i) = value();
      for(size_t j=0;j<A.size2();++j)
        A(i,j) = value();
    }
    ublas::noalias(tst) = ublas::prod(A,x);  
    double tol = 0.01;
    for(size_t i=0;i<A.size1();++i)
    {
      y(i) = big::prod_row( A,x,i);
      BOOST_CHECK_CLOSE( double( y(i) ), double( tst(i) ), tol );
    }
  }
}

BOOST_AUTO_TEST_SUITE_END();
