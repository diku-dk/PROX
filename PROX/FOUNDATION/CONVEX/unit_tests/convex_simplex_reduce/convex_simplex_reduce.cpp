#include <eigenhelperall.h>
#include <convex.h>

#define BOOST_AUTO_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include <boost/test/unit_test_suite.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/test/test_tools.hpp>

BOOST_AUTO_TEST_SUITE(convex_reduce_simplex);

BOOST_AUTO_TEST_CASE(case_by_case_test)
{
    using T = double;
    typedef convex::Simplex<T> simplex_type;

    EigenVector3<T> const zero = EigenVector3<T>(0.0, 0.0, 0.0);

  // Simplex is a single point
    {
        EigenVector3<T> const v1 = EigenVector3<T>(1.0, 1.0, 1.0);

        simplex_type S;

        convex::add_point_to_simplex(v1, v1, v1, S);

        EigenVector3<T> v;
        EigenVector3<T> a;
        EigenVector3<T> b;

        v = convex::reduce_simplex(S, a, b);

        BOOST_CHECK(convex::dimension(S) == 1u);

        BOOST_CHECK(v == zero);
        BOOST_CHECK(a == v1);
        BOOST_CHECK(b == v1);
    }
  // Simplex is an edge that can not be reduced
    {
        EigenVector3<T> const v1 = EigenVector3<T>(1.0, 0.0, 1.0);
        EigenVector3<T> const v2 = EigenVector3<T>(1.0, 0.0, -1.0);

        simplex_type S;

        convex::add_point_to_simplex(v1, v1, v1, S);
        convex::add_point_to_simplex(v2, v2, v2, S);

        EigenVector3<T> v;
        EigenVector3<T> a;
        EigenVector3<T> b;

        v = convex::reduce_simplex(S, a, b);

        BOOST_CHECK(convex::dimension(S) == 2u);

        BOOST_CHECK(v == zero);

        EigenVector3<T> const tst = EigenVector3<T>(1.0, 0.0, 0.0);
        BOOST_CHECK(a == tst);
        BOOST_CHECK(b == tst);
    }
  // Simplex is an edge that can be reduced to a single vertex
    {
        EigenVector3<T> const v1 = EigenVector3<T>(1.0, 1.0, 1.0);
        EigenVector3<T> const v2 = EigenVector3<T>(2.0, 2.0, 2.0);

        simplex_type S;

        convex::add_point_to_simplex(v1, v1, v1, S);
        convex::add_point_to_simplex(v2, v2, v2, S);

        EigenVector3<T> v;
        EigenVector3<T> a;
        EigenVector3<T> b;

        v = convex::reduce_simplex(S, a, b);

        BOOST_CHECK(convex::dimension(S) == 1u);

        BOOST_CHECK(v == zero);
        BOOST_CHECK(a == v1);
        BOOST_CHECK(b == v1);
    }
  // Simplex is a triangle that can not be reduced
    {
        EigenVector3<T> const v1 = EigenVector3<T>(1.0, 0.0, 0.0);
        EigenVector3<T> const v2 = EigenVector3<T>(0.0, 1.0, 0.0);
        EigenVector3<T> const v3 = EigenVector3<T>(0.0, 0.0, 1.0);

        simplex_type S;

        convex::add_point_to_simplex(v1, v1, v1, S);
        convex::add_point_to_simplex(v2, v2, v2, S);
        convex::add_point_to_simplex(v3, v3, v3, S);

        EigenVector3<T> v;
        EigenVector3<T> a;
        EigenVector3<T> b;

        v = convex::reduce_simplex(S, a, b);

        BOOST_CHECK(convex::dimension(S) == 3u);

        BOOST_CHECK(v == zero);
        BOOST_CHECK_CLOSE(a(0), 0.33333333333333331, 0.01);
        BOOST_CHECK_CLOSE(a(1), 0.33333333333333331, 0.01);
        BOOST_CHECK_CLOSE(a(2), 0.33333333333333331, 0.01);
        BOOST_CHECK_CLOSE(b(0), 0.33333333333333331, 0.01);
        BOOST_CHECK_CLOSE(b(1), 0.33333333333333331, 0.01);
        BOOST_CHECK_CLOSE(b(2), 0.33333333333333331, 0.01);
    }
  // Simplex is a triangle that can be reduced to an edge
    {
        EigenVector3<T> const v1 = EigenVector3<T>(1.0, -1.0, 0.0);
        EigenVector3<T> const v2 = EigenVector3<T>(2.0, 0.0, 0.0);
        EigenVector3<T> const v3 = EigenVector3<T>(1.0, 1.0, 0.0);

        simplex_type S;

        convex::add_point_to_simplex(v1, v1, v1, S);
        convex::add_point_to_simplex(v2, v2, v2, S);
        convex::add_point_to_simplex(v3, v3, v3, S);

        EigenVector3<T> v;
        EigenVector3<T> a;
        EigenVector3<T> b;

        v = convex::reduce_simplex(S, a, b);

        BOOST_CHECK(convex::dimension(S) == 2u);

        BOOST_CHECK(v == zero);
        BOOST_CHECK_CLOSE(a(0), 1.0, 0.01);
        BOOST_CHECK_CLOSE(a(1), 0.0, 0.01);
        BOOST_CHECK_CLOSE(a(2), 0.0, 0.01);
        BOOST_CHECK_CLOSE(b(0), 1.0, 0.01);
        BOOST_CHECK_CLOSE(b(1), 0.0, 0.01);
        BOOST_CHECK_CLOSE(b(2), 0.0, 0.01);
    }
  // Simplex is a triangle that can be reduced to a point
    {
        EigenVector3<T> const v1 = EigenVector3<T>(2.0, -1.0, 0.0);
        EigenVector3<T> const v2 = EigenVector3<T>(2.0, 1.0, 0.0);
        EigenVector3<T> const v3 = EigenVector3<T>(1.0, 0.0, 0.0);

        simplex_type S;

        convex::add_point_to_simplex(v1, v1, v1, S);
        convex::add_point_to_simplex(v2, v2, v2, S);
        convex::add_point_to_simplex(v3, v3, v3, S);

        EigenVector3<T> v;
        EigenVector3<T> a;
        EigenVector3<T> b;

        v = convex::reduce_simplex(S, a, b);

        BOOST_CHECK(convex::dimension(S) == 1u);

        BOOST_CHECK(v == zero);
        BOOST_CHECK_CLOSE(a(0), 1.0, 0.01);
        BOOST_CHECK_CLOSE(a(1), 0.0, 0.01);
        BOOST_CHECK_CLOSE(a(2), 0.0, 0.01);
        BOOST_CHECK_CLOSE(b(0), 1.0, 0.01);
        BOOST_CHECK_CLOSE(b(1), 0.0, 0.01);
        BOOST_CHECK_CLOSE(b(2), 0.0, 0.01);
    }
  // Origin is inside a tetrahedron, so tetrahedron can not be reduced
    {
        EigenVector3<T> const v1 = EigenVector3<T>(-1.0, 0.0, -1.0);
        EigenVector3<T> const v2 = EigenVector3<T>(1.0, 0.0, -1.0);
        EigenVector3<T> const v3 = EigenVector3<T>(0.0, -1.0, 1.0);
        EigenVector3<T> const v4 = EigenVector3<T>(0.0, 1.0, 1.0);

        simplex_type S;

        convex::add_point_to_simplex(v1, v1, v1, S);
        convex::add_point_to_simplex(v2, v2, v2, S);
        convex::add_point_to_simplex(v3, v3, v3, S);
        convex::add_point_to_simplex(v4, v4, v4, S);

        EigenVector3<T> v;
        EigenVector3<T> a;
        EigenVector3<T> b;

        v = convex::reduce_simplex(S, a, b);

        BOOST_CHECK(convex::dimension(S) == 4u);

        BOOST_CHECK(v == zero);
        BOOST_CHECK(a == v);
        BOOST_CHECK(b == v);
    }
  // Simplex is a tetrahedron that can be reduced to a triangle
    {
        EigenVector3<T> const v1 = EigenVector3<T>(1.0, 0.0, 0.0);
        EigenVector3<T> const v2 = EigenVector3<T>(0.0, 1.0, 0.0);
        EigenVector3<T> const v3 = EigenVector3<T>(0.0, 0.0, 1.0);
        EigenVector3<T> const v4 = EigenVector3<T>(2.0, 2.0, 2.0);

        simplex_type S;

        convex::add_point_to_simplex(v1, v1, v1, S);
        convex::add_point_to_simplex(v2, v2, v2, S);
        convex::add_point_to_simplex(v3, v3, v3, S);
        convex::add_point_to_simplex(v4, v4, v4, S);

        EigenVector3<T> v;
        EigenVector3<T> a;
        EigenVector3<T> b;

        v = convex::reduce_simplex(S, a, b);

        BOOST_CHECK(convex::dimension(S) == 3u);

        BOOST_CHECK(v == zero);
        BOOST_CHECK_CLOSE(a(0), 0.33333333333333331, 0.01);
        BOOST_CHECK_CLOSE(a(1), 0.33333333333333331, 0.01);
        BOOST_CHECK_CLOSE(a(2), 0.33333333333333331, 0.01);
        BOOST_CHECK_CLOSE(b(0), 0.33333333333333331, 0.01);
        BOOST_CHECK_CLOSE(b(1), 0.33333333333333331, 0.01);
        BOOST_CHECK_CLOSE(b(2), 0.33333333333333331, 0.01);
    }
  // Simplex is a tetrahedron that can be reduced to a point
    {
        EigenVector3<T> const v1 = EigenVector3<T>(1.0, 0.0, 0.0);
        EigenVector3<T> const v2 = EigenVector3<T>(0.0, 1.0, 0.0);
        EigenVector3<T> const v3 = EigenVector3<T>(0.0, 0.0, 1.0);
        EigenVector3<T> const v4 = EigenVector3<T>(0.25, 0.25, 0.25);

        simplex_type S;

        convex::add_point_to_simplex(v1, v1, v1, S);
        convex::add_point_to_simplex(v2, v2, v2, S);
        convex::add_point_to_simplex(v3, v3, v3, S);
        convex::add_point_to_simplex(v4, v4, v4, S);

        EigenVector3<T> v;
        EigenVector3<T> a;
        EigenVector3<T> b;

        v = convex::reduce_simplex(S, a, b);

        BOOST_CHECK(convex::dimension(S) == 1u);

        BOOST_CHECK(v == zero);
        BOOST_CHECK(a == v4);
        BOOST_CHECK(b == v4);
    }
  // Simplex is a tetrahedron that can be reduced to an edge
    {
        EigenVector3<T> const v1 = EigenVector3<T>(1.0, -1.0, 0.0);
        EigenVector3<T> const v2 = EigenVector3<T>(1.0, 1.0, 0.0);
        EigenVector3<T> const v3 = EigenVector3<T>(2.0, 0.0, -1.0);
        EigenVector3<T> const v4 = EigenVector3<T>(2.0, 0.0, 1.0);

        simplex_type S;

        convex::add_point_to_simplex(v1, v1, v1, S);
        convex::add_point_to_simplex(v2, v2, v2, S);
        convex::add_point_to_simplex(v3, v3, v3, S);
        convex::add_point_to_simplex(v4, v4, v4, S);

        EigenVector3<T> v;
        EigenVector3<T> a;
        EigenVector3<T> b;

        v = convex::reduce_simplex(S, a, b);

        BOOST_CHECK(convex::dimension(S) == 2u);

        BOOST_CHECK(v == zero);
        BOOST_CHECK_CLOSE(a(0), 1.0, 0.01);
        BOOST_CHECK_CLOSE(a(1), 0.0, 0.01);
        BOOST_CHECK_CLOSE(a(2), 0.0, 0.01);
        BOOST_CHECK_CLOSE(b(0), 1.0, 0.01);
        BOOST_CHECK_CLOSE(b(1), 0.0, 0.01);
        BOOST_CHECK_CLOSE(b(2), 0.0, 0.01);
    }
}

BOOST_AUTO_TEST_SUITE_END();
