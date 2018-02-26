#include <big_forward_gauss_seidel.h>
#include <boost/cast.hpp>

#define BOOST_AUTO_TEST_MAIN
#include <boost/test/auto_unit_test.hpp>
#include <boost/test/unit_test_suite.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/test_tools.hpp>

template<typename matrix_type,typename vector_type>
void test(matrix_type const & A, vector_type  & x, vector_type const & b, vector_type const & y)
{
  typedef typename matrix_type::value_type real_type;
  typedef typename matrix_type::size_type  size_type;

  real_type const tol = boost::numeric_cast<real_type>(0.01);

  big::forward_gauss_seidel(A,x,b);

  for(size_type i = 0; i < x.size();++i)
    BOOST_CHECK_CLOSE( real_type( x(i) ), real_type( y(i) ), tol );
}


BOOST_AUTO_TEST_SUITE(big_forward_gauss_seidel);

BOOST_AUTO_TEST_CASE(logic_and_valid_arguments_testing)
{
  typedef ublas::compressed_matrix<double> matrix_type;
  typedef ublas::vector<double>            vector_type;

  {
    matrix_type A;
    vector_type x,b;
    BOOST_CHECK_THROW( big::forward_gauss_seidel(A,x,b), std::invalid_argument );
  }
  {
    matrix_type A;
    vector_type x,b;
    A.resize(4,4,false);
    BOOST_CHECK_THROW( big::forward_gauss_seidel(A,x,b), std::invalid_argument );
  }
  {
    matrix_type A;
    vector_type x,b;
    A.resize(4,4,false);
    b.resize(4,false);
    BOOST_CHECK_THROW( big::forward_gauss_seidel(A,x,b), std::invalid_argument );
  }
  {
    matrix_type A;
    vector_type x,b;
    A.resize(4,4,false);
    x.resize(4,false);
    BOOST_CHECK_THROW( big::forward_gauss_seidel(A,x,b), std::invalid_argument );
  }
  {
    matrix_type A;
    vector_type x,b;
    A.resize(4,4,false);
    x.resize(4,false);
    b.resize(4,false);
    A(0,0) = A(1,1) = A(2,2) = A(3,3) = 1.0;
    BOOST_CHECK_NO_THROW( big::forward_gauss_seidel(A,x,b) );
  }
  {
    matrix_type A;
    vector_type x,b;
    x.resize(4,false);
    b.resize(4,false);
    BOOST_CHECK_THROW( big::forward_gauss_seidel(A,x,b), std::invalid_argument );
  }
  {
    matrix_type A;
    vector_type x,b;
    A.resize(4,4,false);
    x.resize(3,false);
    b.resize(4,false);
    BOOST_CHECK_THROW( big::forward_gauss_seidel(A,x,b), std::invalid_argument );
  }
  {
    matrix_type A;
    vector_type x,b;
    A.resize(4,4,false);
    x.resize(4,false);
    b.resize(3,false);
    BOOST_CHECK_THROW( big::forward_gauss_seidel(A,x,b), std::invalid_argument );
  }
  {
    matrix_type A;
    vector_type x,b;
    A.resize(4,4,false);
    x.resize(3,false);
    b.resize(3,false);
    BOOST_CHECK_THROW( big::forward_gauss_seidel(A,x,b), std::invalid_argument );
  }
}


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
  y(0) = -25.000000;
  y(1) = 12.500000;
  y(2) = 4.909091;
  y(3) = 4.835227;
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
  test(A,x,b,y);

  x(0) = 0.133773;
  x(1) = 0.207133;
  x(2) = 0.607199;
  x(3) = 0.629888;
  b(0) = 0.586918;
  b(1) = 0.057581;
  b(2) = 0.367568;
  b(3) = 0.631451;
  y(0) = -0.585224;
  y(1) = 1.048565;
  y(2) = 17.204089;
  y(3) = -54.090245;
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
  test(A,x,b,y);

  x(0) = 0.717634;
  x(1) = 0.692669;
  x(2) = 0.084079;
  x(3) = 0.454355;
  b(0) = 0.804872;
  b(1) = 0.908398;
  b(2) = 0.231894;
  b(3) = 0.239313;
  y(0) = 1.821684;
  y(1) = 1.248641;
  y(2) = 0.323928;
  y(3) = 1.029964;
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
  test(A,x,b,y);

  x(0) = 0.297406;
  x(1) = 0.049162;
  x(2) = 0.693180;
  x(3) = 0.650106;
  b(0) = 0.193893;
  b(1) = 0.904812;
  b(2) = 0.569206;
  b(3) = 0.631790;
  y(0) = 0.197249;
  y(1) = 1.085136;
  y(2) = 0.103801;
  y(3) = 0.665095;
  A(0,0) = 0.982988;
  A(0,1) = 0.000000;
  A(0,2) = 0.000000;
  A(0,3) = 0.000000;
  A(1,0) = 0.552673;
  A(1,1) = 0.733363;
  A(1,2) = 0.000000;
  A(1,3) = 0.000000;
  A(2,0) = 0.400074;
  A(2,1) = 0.375885;
  A(2,2) = 0.793872;
  A(2,3) = 0.000000;
  A(3,0) = 0.198789;
  A(3,1) = 0.009876;
  A(3,2) = 0.919957;
  A(3,3) = 0.731277;
  test(A,x,b,y);

  x(0) = 0.234413;
  x(1) = 0.548782;
  x(2) = 0.931583;
  x(3) = 0.335197;
  b(0) = 0.185904;
  b(1) = 0.700635;
  b(2) = 0.982709;
  b(3) = 0.806638;
  y(0) = -0.943403;
  y(1) = 0.321753;
  y(2) = 1.214908;
  y(3) = 1.648401;
  A(0,0) = 0.655531;
  A(0,1) = 0.397184;
  A(0,2) = 0.371608;
  A(0,3) = 0.716542;
  A(1,0) = 0.000000;
  A(1,1) = 0.413629;
  A(1,2) = 0.425253;
  A(1,3) = 0.511311;
  A(2,0) = 0.000000;
  A(2,1) = 0.000000;
  A(2,2) = 0.594663;
  A(2,3) = 0.776401;
  A(3,0) = 0.000000;
  A(3,1) = 0.000000;
  A(3,2) = 0.000000;
  A(3,3) = 0.489345;
  test(A,x,b,y);
}

BOOST_AUTO_TEST_SUITE_END();
