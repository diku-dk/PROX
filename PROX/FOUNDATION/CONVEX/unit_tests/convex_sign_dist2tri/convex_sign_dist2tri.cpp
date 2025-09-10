#include <tiny.h>
#include <convex.h>

#define BOOST_AUTO_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include <boost/test/unit_test_suite.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/test/test_tools.hpp>

BOOST_AUTO_TEST_SUITE(convex_signed_distance_to_triangle);

BOOST_AUTO_TEST_CASE(case_by_case_test)
{
    using T = double;
    using real_type = T;

    EigenVector3<T> a = EigenVector3<T>(0.0, 0.0, 0.0);
    EigenVector3<T> b = EigenVector3<T>(1.0, 0.0, 0.0);
    ;
    EigenVector3<T> c = EigenVector3<T>(0.0, 1.0, 0.0);
    ;
    EigenVector3<T> q = EigenVector3<T>(0.33, 0.33, -1.0);

    {
        EigenVector3<T> p = EigenVector3<T>(0.33, 0.33, 1.0);
        real_type sign_p = 0.0;
        sign_p = convex::signed_distance_to_triangle(p, a, b, c, q);
        BOOST_CHECK_CLOSE(sign_p, 1.0, 0.01);
    }

    {
        EigenVector3<T> p = EigenVector3<T>(0.1, 0.1, -1.0);
        real_type sign_p = 0.0;
        sign_p = convex::signed_distance_to_triangle(p, a, b, c, q);
        BOOST_CHECK_CLOSE(sign_p, -1.0, 0.01);
    }

    {
        EigenVector3<T> p = EigenVector3<T>(0.1, 0.1, 0.0);
        real_type sign_p = 10.0;
        sign_p = convex::signed_distance_to_triangle(p, a, b, c, q);
        BOOST_CHECK_CLOSE(sign_p, 0.0, 0.01);
    }
}

BOOST_AUTO_TEST_SUITE_END();
