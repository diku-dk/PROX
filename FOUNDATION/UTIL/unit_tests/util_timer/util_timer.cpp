#include <util_timer.h>

#define BOOST_AUTO_TEST_MAIN
#include <boost/test/auto_unit_test.hpp>
#include <boost/test/unit_test_suite.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/test_tools.hpp>

BOOST_AUTO_TEST_SUITE(util);

BOOST_AUTO_TEST_CASE(timer_test)
{
  
  util::Timer timer;
  timer.start();

  BOOST_CHECK( timer() == 0.0f );

  {
    float fake = 1.0f;
    for(size_t i=0u;i< 1000000u;++i)
      fake *= 2.0f;
  }
  
  timer.stop();

  BOOST_CHECK( timer() > 0.0f );  
}

BOOST_AUTO_TEST_SUITE_END();
