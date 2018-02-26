#include <util_file_locator.h>

#define BOOST_AUTO_TEST_MAIN
#include <boost/test/auto_unit_test.hpp>
#include <boost/test/unit_test_suite.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/test_tools.hpp>

BOOST_AUTO_TEST_SUITE(util);

BOOST_AUTO_TEST_CASE(file_locator_test)
{
  util::FileLocator myfile1("/Users/kenny/Documents/PROX_TRUNK/resources/objs", "bunny.obj");
  util::FileLocator myfile2("/Users/kenny/Documents/PROX_TRUNK/resources/objs/", "bunny.obj");
  util::FileLocator myfile3("/Users/kenny/Documents/PROX_TRUNK/resources/objs/", "/bunny.obj");
  util::FileLocator myfile4("/Users/kenny/Documents/PROX_TRUNK/resources/objs/bunny.obj");
  util::FileLocator myfile5("", "bunny.obj");
  util::FileLocator myfile6("/Users/kenny/Documents/PROX_TRUNK/resources/textures", "bunny.obj");
  util::FileLocator myfile7("/Users/kenny/Documents/PROX_TRUNK/resources/", "objs/");

  BOOST_CHECK_NO_THROW( myfile1.get_filename() );
  BOOST_CHECK_NO_THROW( myfile2.get_filename() );
  BOOST_CHECK_NO_THROW( myfile3.get_filename() );
  BOOST_CHECK_NO_THROW( myfile4.get_filename() );

  BOOST_CHECK_THROW(myfile5.get_filename(), std::logic_error );
  BOOST_CHECK_THROW(myfile6.get_filename(), std::logic_error );
  BOOST_CHECK_THROW(myfile7.get_filename(), std::logic_error );

  BOOST_CHECK_NO_THROW( util::get_full_filename("/Users/kenny/Documents/PROX_TRUNK/resources/objs", "bunny.obj") );
  BOOST_CHECK_NO_THROW( util::get_full_filename("/Users/kenny/Documents/PROX_TRUNK/resources/objs/", "bunny.obj") );
  BOOST_CHECK_NO_THROW( util::get_full_filename("/Users/kenny/Documents/PROX_TRUNK/resources/objs/", "/bunny.obj") );
  BOOST_CHECK_NO_THROW( util::get_full_filename("", "/Users/kenny/Documents/PROX_TRUNK/resources/objs/bunny.obj") );

  BOOST_CHECK_THROW( util::get_full_filename("", "bunny.obj"), std::logic_error );
  BOOST_CHECK_THROW( util::get_full_filename("/Users/kenny/Documents/PROX_TRUNK/resources/textures", "bunny.obj"), std::logic_error );
  BOOST_CHECK_THROW( util::get_full_filename("/Users/kenny/Documents/PROX_TRUNK/resources/", "objs/"), std::logic_error );
}

BOOST_AUTO_TEST_SUITE_END();
