#include <tiny.h>
#include <convex.h>

#define BOOST_AUTO_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include <boost/test/unit_test_suite.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/test/test_tools.hpp>

BOOST_AUTO_TEST_SUITE(convex_outside_triangle);

BOOST_AUTO_TEST_CASE(case_by_case_test)
{
    using T = double;

    const EigenVector3<T> a = EigenVector3<T>(0.0, 0.0, 0.0);
    const EigenVector3<T> b = EigenVector3<T>(1.0, 0.0, 0.0);
    ;
    const EigenVector3<T> c = EigenVector3<T>(0.0, 1.0, 0.0);
    ;
    const EigenVector3<T> q = EigenVector3<T>(0.33, 0.33, -1.0);

    {
        EigenVector3<T> p = EigenVector3<T>(0.33, 0.33, 1.0);
        bool outside = convex::outside_triangle(p, a, b, c, q);
        BOOST_CHECK(outside);
    }
    {
        EigenVector3<T> p = EigenVector3<T>(0.1, 0.1, -1.0);
        bool outside = convex::outside_triangle(p, a, b, c, q);
        BOOST_CHECK(!outside);
    }
    {
        EigenVector3<T> p = EigenVector3<T>(0.1, 0.1, 0.0);
        bool outside = convex::outside_triangle(p, a, b, c, q);
        BOOST_CHECK(outside);
    }
}

BOOST_AUTO_TEST_SUITE_END();
