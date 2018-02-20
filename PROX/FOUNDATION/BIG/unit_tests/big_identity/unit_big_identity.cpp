#include <big_types.h>
#include <big_identity.h>
#include <boost/cast.hpp>

#define BOOST_AUTO_TEST_MAIN
#include <boost/test/auto_unit_test.hpp>
#include <boost/test/unit_test_suite.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/test_tools.hpp>

template<typename matrix_type,typename vector_type>
void test(matrix_type const & A, vector_type  & x, vector_type const & b)
{
  typedef typename matrix_type::value_type real_type;
  typedef typename matrix_type::size_type  size_type;

  real_type const tol = boost::numeric_cast<real_type>(0.01);

  big::identity(A,x,b);   // solves A x = b assuming A = I eventhough it is not!
  for(size_type i = 0; i < x.size();++i)
    BOOST_CHECK_CLOSE( real_type( x(i) ), real_type( b(i) ), tol );
  
  x.clear();
  
  big::identity(A,x,b);
  for(size_type i = 0; i < x.size();++i)
    BOOST_CHECK_CLOSE( real_type( x(i) ), real_type( b(i) ), tol );
  
}


BOOST_AUTO_TEST_SUITE(big_identity);

BOOST_AUTO_TEST_CASE(correctness_testing)
{
  typedef ublas::compressed_matrix<double> matrix_type;
  typedef ublas::vector<double>            vector_type;

  matrix_type A;
  A.resize(4,4,false);
  vector_type x,b,y;
  x.resize(4,false);
  b.resize(4,false);
  y.resize(4,false);

  x(0) = 1.000000;
  x(1) = 2.000000;
  x(2) = 3.000000;
  x(3) = 4.000000;
  
  b(0) = 4.000000;
  b(1) = 3.000000;
  b(2) = 2.000000;
  b(3) = 1.000000;
  
  A(0,0) = 1.000000;
  A(0,1) = 2.000000;
  A(0,2) = 3.000000;
  A(0,3) = 4.000000;
  A(1,0) = 5.000000;
  A(1,1) = 6.000000;
  A(1,2) = 7.000000;
  A(1,3) = 8.000000;
  A(2,0) = 9.000000;
  A(2,1) = 10.000000;
  A(2,2) = 11.000000;
  A(2,3) = 12.000000;
  A(3,0) = 13.000000;
  A(3,1) = 14.000000;
  A(3,2) = 15.000000;
  A(3,3) = 16.000000;
  test(A,x,b);

  x(0) = 0.133773;
  x(1) = 0.207133;
  x(2) = 0.607199;
  x(3) = 0.629888;
  
  b(0) = 0.586918;
  b(1) = 0.057581;
  b(2) = 0.367568;
  b(3) = 0.631451;
  
  A(0,0) = 0.370477;
  A(0,1) = 0.027185;
  A(0,2) = 0.683116;
  A(0,3) = 0.608540;
  A(1,0) = 0.575148;
  A(1,1) = 0.312685;
  A(1,2) = 0.092842;
  A(1,3) = 0.015760;
  A(2,0) = 0.451425;
  A(2,1) = 0.012863;
  A(2,2) = 0.035338;
  A(2,3) = 0.016355;
  A(3,0) = 0.043895;
  A(3,1) = 0.383967;
  A(3,2) = 0.612395;
  A(3,3) = 0.190075;
  test(A,x,b);

  x(0) = 0.717634;
  x(1) = 0.692669;
  x(2) = 0.084079;
  x(3) = 0.454355;
  
  b(0) = 0.804872;
  b(1) = 0.908398;
  b(2) = 0.231894;
  b(3) = 0.239313;
  
  A(0,0) = 0.441828;
  A(0,1) = 0.000000;
  A(0,2) = 0.000000;
  A(0,3) = 0.000000;
  A(1,0) = 0.000000;
  A(1,1) = 0.727509;
  A(1,2) = 0.000000;
  A(1,3) = 0.000000;
  A(2,0) = 0.000000;
  A(2,1) = 0.000000;
  A(2,2) = 0.715883;
  A(2,3) = 0.000000;
  A(3,0) = 0.000000;
  A(3,1) = 0.000000;
  A(3,2) = 0.000000;
  A(3,3) = 0.232350;
  test(A,x,b);
}

BOOST_AUTO_TEST_SUITE_END();
