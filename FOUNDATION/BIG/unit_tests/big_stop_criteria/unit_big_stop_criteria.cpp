#include <big_types.h>
#include <big_absolute_convergence.h>
#include <big_relative_convergence.h>
#include <big_stagnation.h>

#define BOOST_AUTO_TEST_MAIN
#include <boost/test/auto_unit_test.hpp>
#include <boost/test/unit_test_suite.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/test_tools.hpp>

BOOST_AUTO_TEST_SUITE(big_stop_criteria);

BOOST_AUTO_TEST_CASE(absolute_testing)
{
  BOOST_CHECK(  big::absolute_convergence(  1.0,   0.01) == false  );
  BOOST_CHECK(  big::absolute_convergence(  0.01,  0.01) == true   );
  BOOST_CHECK(  big::absolute_convergence(  0.005, 0.01) == true   );
  BOOST_CHECK(  big::absolute_convergence(  0.0,   0.01) == true   );
  BOOST_CHECK(  big::absolute_convergence(  0.0,   0.0 ) == true   );

//  BOOST_CHECK_THROW(  big::absolute_convergence( -1.0,   0.01), std::invalid_argument  ); // 2012-09-30 Kenny: Should we throw an exception?
//  BOOST_CHECK_THROW(  big::absolute_convergence(  1.0,  -0.01), std::invalid_argument  ); // 2012-09-30 Kenny: Should we throw an exception?
  
}

BOOST_AUTO_TEST_CASE(relative_testing)
{
  
  // Test is:  fabs(old - new) / fabs(old) <= tol
  
  BOOST_CHECK(  big::relative_convergence(  1.0, 1.0, 0.01 ) == true   );
  BOOST_CHECK(  big::relative_convergence(  1.0, 1.0, 0.0  ) == true   );
  
  BOOST_CHECK(  big::relative_convergence(  1.0, 1.1, 0.11 ) == true   );
  BOOST_CHECK(  big::relative_convergence(  1.0, 1.1, 0.09 ) == false  );
  BOOST_CHECK(  big::relative_convergence(  1.0, 0.9, 0.11 ) == true   );
  BOOST_CHECK(  big::relative_convergence(  1.0, 0.9, 0.09 ) == false  );

  BOOST_CHECK(  big::relative_convergence(  0.0, 0.0, 0.11 ) == true   );
  BOOST_CHECK(  big::relative_convergence(  0.0, .1, 0.11  ) == false  );
  
  
  //  BOOST_CHECK_THROW(  big::relative_convergence( 1.0, 1.0, -0.01), std::invalid_argument  ); // 2012-09-30 Kenny: Should we throw an exception?
}

BOOST_AUTO_TEST_CASE(stagnation_testing)
{
  typedef ublas::vector<double>            vector_type;
  vector_type x;
  vector_type y;
  
  
  x.resize(2);
  y.resize(2);
  
  x(0) = 0.01;
  x(1) = 0.01;
  
  y(0) = 0.01;
  y(1) = 0.01;

  BOOST_CHECK(  big::stagnation(  x, y, 0.01 ) == true   );

  
  x(0) = 0.01;
  x(1) = 0.01;
  
  y(0) = 0.01 + 2*0.01;
  y(1) = 0.01;
  
  BOOST_CHECK(  big::stagnation(  x, y, 0.01 ) == false   );
  
  BOOST_CHECK(  big::stagnation(  x, y, 0.0  ) == false  );
  BOOST_CHECK(  big::stagnation(  x, x, 0.0  ) == true   );

  //  BOOST_CHECK_THROW(  big::stagnation( x, y, -0.01), std::invalid_argument  ); // 2012-09-30 Kenny: Should we throw an exception?

}

BOOST_AUTO_TEST_SUITE_END();
