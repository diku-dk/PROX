#include <big_types.h>
#include <big_project.h>

#define BOOST_AUTO_TEST_MAIN
#include <boost/test/auto_unit_test.hpp>
#include <boost/test/unit_test_suite.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/test_tools.hpp>

BOOST_AUTO_TEST_SUITE(big_project);

BOOST_AUTO_TEST_CASE(test_case)
{
  ublas::vector<double> x;
  
  x.resize(10,false);
  
  x(0) =  2.0;
  x(1) = -3.0;
  x(2) =  4.0;
  x(3) = -5.0;
  x(4) =  6.0;
  x(5) = -7.0;
  x(6) =  8.0;
  x(7) = -9.0;
  x(8) =  10.0;
  x(9) = -11.0;
  
  double const tol =0.01;
  ublas::vector<double>  y = x;
  
  big::project(x,y);
  
  BOOST_CHECK_CLOSE( y(0) , 2.0 , tol );
  BOOST_CHECK_CLOSE( y(1) , 0.0 , tol );
  BOOST_CHECK_CLOSE( y(2) , 4.0 , tol );
  BOOST_CHECK_CLOSE( y(3) , 0.0 , tol );
  BOOST_CHECK_CLOSE( y(4) , 6.0 , tol );
  BOOST_CHECK_CLOSE( y(5) , 0.0 , tol );
  BOOST_CHECK_CLOSE( y(6) , 8.0 , tol );
  BOOST_CHECK_CLOSE( y(7) , 0.0 , tol );
  BOOST_CHECK_CLOSE( y(8) , 10.0 , tol );
  BOOST_CHECK_CLOSE( y(9) , 0.0 , tol );
  
}

BOOST_AUTO_TEST_SUITE_END();
