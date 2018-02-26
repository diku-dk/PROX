#include <convex.h>

#define BOOST_AUTO_TEST_MAIN
#include <boost/test/auto_unit_test.hpp>
#include <boost/test/unit_test_suite.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/test_tools.hpp>

#include <string>

BOOST_AUTO_TEST_SUITE(convex_constants);

BOOST_AUTO_TEST_CASE(simple_testing)
{
  std::string msg = "";

  BOOST_CHECK_NO_THROW( msg = convex::get_status_message( convex::ABSOLUTE_CONVERGENCE ) );
  BOOST_CHECK( msg == "Absolute convergence test passed" );

  BOOST_CHECK_NO_THROW( msg = convex::get_status_message( convex::EXCEEDED_MAX_ITERATIONS_LIMIT ) );
  BOOST_CHECK( msg == "Maximum iteration limit was exceeded" );

  BOOST_CHECK_NO_THROW( msg = convex::get_status_message( convex::INTERSECTION ) );
  BOOST_CHECK( msg == "Intersection was found" );

  BOOST_CHECK_NO_THROW( msg = convex::get_status_message( convex::ITERATING ) );
  BOOST_CHECK( msg == "Unexpected termination while iterating" );

  BOOST_CHECK_NO_THROW( msg = convex::get_status_message( convex::NON_DESCEND_DIRECTION ) );
  BOOST_CHECK( msg == "Non descent direction was encountered" );

  BOOST_CHECK_NO_THROW( msg = convex::get_status_message( convex::RELATIVE_CONVERGENCE ) );
  BOOST_CHECK( msg == "Relative convergence test passed" );

  BOOST_CHECK_NO_THROW( msg = convex::get_status_message( convex::STAGNATION ) );
  BOOST_CHECK( msg == "Stagnation test passed" );

  BOOST_CHECK_NO_THROW( msg = convex::get_status_message( convex::SIMPLEX_EXPANSION_FAILED ) );
  BOOST_CHECK( msg == "Simplex expansion failure" );

  BOOST_CHECK_NO_THROW( msg = convex::get_status_message( 7654u ) );
  BOOST_CHECK( msg == "unrecognised error" );
}

BOOST_AUTO_TEST_SUITE_END();
