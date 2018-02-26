#include <util_folder_locator.h>

#define BOOST_AUTO_TEST_MAIN
#include <boost/test/auto_unit_test.hpp>
#include <boost/test/unit_test_suite.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/test_tools.hpp>

BOOST_AUTO_TEST_SUITE(util);

BOOST_AUTO_TEST_CASE(file_locator_test)
{
  util::FolderLocator f1("/Users/kenny/Documents/PROX_TRUNK/resources/", "objs");
  util::FolderLocator f2("/Users/kenny/Documents/PROX_TRUNK/resources", "objs");
  util::FolderLocator f3("/Users/kenny/Documents/PROX_TRUNK/resources/", "/objs");
  util::FolderLocator f4("/Users/kenny/Documents/PROX_TRUNK/resources/objs/");
  util::FolderLocator f5("/Users/kenny/Documents/PROX_TRUNK/resources/objs");

  util::FolderLocator f6("", "objs");
  util::FolderLocator f7("/Users/kenny/Documents/PROX_TRUNK/bin", "objs");
  util::FolderLocator f8("/Users/kenny/Documents/PROX_TRUNK/resources/objs/", "bunny.obj");


  BOOST_CHECK_NO_THROW( f1.get_folder_name() );
  BOOST_CHECK_NO_THROW( f2.get_folder_name() );
  BOOST_CHECK_NO_THROW( f3.get_folder_name() );
  BOOST_CHECK_NO_THROW( f4.get_folder_name() );
  BOOST_CHECK_NO_THROW( f5.get_folder_name() );

  BOOST_CHECK_THROW(f6.get_folder_name(), std::logic_error );
  BOOST_CHECK_THROW(f7.get_folder_name(), std::logic_error );
  BOOST_CHECK_THROW(f8.get_folder_name(), std::logic_error );

}

BOOST_AUTO_TEST_SUITE_END();
