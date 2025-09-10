#include <geometry.h>
#include <eigenhelperall.h>

#define BOOST_AUTO_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include <boost/test/unit_test_suite.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/test/test_tools.hpp>

#include <vector>

BOOST_AUTO_TEST_SUITE(geometry);

BOOST_AUTO_TEST_CASE(make_dops)
{
    using T = float;

    {
        geometry::DOP<T, 6> A;
        BOOST_CHECK(geometry::is_valid(A) == false);
    }

    {
        geometry::DirectionTable<T, 3> DT = geometry::make3<T>();

        std::vector<EigenVector3<T>> points(8);

        points[0] = EigenVector3<T>(-1.0f, -1.0f, -1.0f);
        points[1] = EigenVector3<T>(1.0f, -1.0f, -1.0f);
        points[2] = EigenVector3<T>(1.0f, 1.0f, -1.0f);
        points[3] = EigenVector3<T>(-1.0f, 1.0f, -1.0f);
        points[4] = EigenVector3<T>(-1.0f, -1.0f, 1.0f);
        points[5] = EigenVector3<T>(1.0f, -1.0f, 1.0f);
        points[6] = EigenVector3<T>(1.0f, 1.0f, 1.0f);
        points[7] = EigenVector3<T>(-1.0f, 1.0f, 1.0f);

        geometry::DOP<T, 6> A = geometry::make_dop(&points[0], &points[0] + 8, DT);

        BOOST_CHECK_EQUAL(A.size(), 6u);

        BOOST_CHECK_EQUAL(A(0).lower(), -1.0);
        BOOST_CHECK_EQUAL(A(1).lower(), -1.0);
        BOOST_CHECK_EQUAL(A(2).lower(), -1.0);
        BOOST_CHECK_EQUAL(A(0).upper(), 1.0);
        BOOST_CHECK_EQUAL(A(1).upper(), 1.0);
        BOOST_CHECK_EQUAL(A(2).upper(), 1.0);
    }

  // Separated kDOPs
    {
        geometry::DirectionTable<T, 3> DT = geometry::make3<T>();

        std::vector<EigenVector3<T>> points(8);

        points[0] = EigenVector3<T>(-1.0f, -1.0f, -1.0f)
                  + EigenVector3<T>(-2.0f, 0.0f, 0.0f);
        points[1] = EigenVector3<T>(1.0f, -1.0f, -1.0f)
                  + EigenVector3<T>(-2.0f, 0.0f, 0.0f);
        points[2] = EigenVector3<T>(1.0f, 1.0f, -1.0f)
                  + EigenVector3<T>(-2.0f, 0.0f, 0.0f);
        points[3] = EigenVector3<T>(-1.0f, 1.0f, -1.0f)
                  + EigenVector3<T>(-2.0f, 0.0f, 0.0f);
        points[4] = EigenVector3<T>(-1.0f, -1.0f, 1.0f)
                  + EigenVector3<T>(-2.0f, 0.0f, 0.0f);
        points[5] = EigenVector3<T>(1.0f, -1.0f, 1.0f)
                  + EigenVector3<T>(-2.0f, 0.0f, 0.0f);
        points[6] = EigenVector3<T>(1.0f, 1.0f, 1.0f)
                  + EigenVector3<T>(-2.0f, 0.0f, 0.0f);
        points[7] = EigenVector3<T>(-1.0f, 1.0f, 1.0f)
                  + EigenVector3<T>(-2.0f, 0.0f, 0.0f);

        geometry::DOP<T, 6> A = geometry::make_dop(&points[0], &points[0] + 8, DT);

        points[0] = EigenVector3<T>(-1.0f, -1.0f, -1.0f)
                  + EigenVector3<T>(2.0f, 0.0f, 0.0f);
        points[1] = EigenVector3<T>(1.0f, -1.0f, -1.0f)
                  + EigenVector3<T>(2.0f, 0.0f, 0.0f);
        points[2] = EigenVector3<T>(1.0f, 1.0f, -1.0f)
                  + EigenVector3<T>(2.0f, 0.0f, 0.0f);
        points[3] = EigenVector3<T>(-1.0f, 1.0f, -1.0f)
                  + EigenVector3<T>(2.0f, 0.0f, 0.0f);
        points[4] = EigenVector3<T>(-1.0f, -1.0f, 1.0f)
                  + EigenVector3<T>(2.0f, 0.0f, 0.0f);
        points[5] = EigenVector3<T>(1.0f, -1.0f, 1.0f)
                  + EigenVector3<T>(2.0f, 0.0f, 0.0f);
        points[6] = EigenVector3<T>(1.0f, 1.0f, 1.0f)
                  + EigenVector3<T>(2.0f, 0.0f, 0.0f);
        points[7] = EigenVector3<T>(-1.0f, 1.0f, 1.0f)
                  + EigenVector3<T>(2.0f, 0.0f, 0.0f);

        geometry::DOP<T, 6> B = geometry::make_dop(&points[0], &points[0] + 8, DT);

        BOOST_CHECK(geometry::overlap_dop_dop(A, B) == false);
        BOOST_CHECK(geometry::overlap_dop_dop(B, A) == false);
    }
  // Penetrating kDOPs
    {
        geometry::DirectionTable<T, 3> DT = geometry::make3<T>();

        std::vector<EigenVector3<T>> points(8);

        points[0] = EigenVector3<T>(-1.0f, -1.0f, -1.0f)
                  + EigenVector3<T>(-0.5f, 0.0f, 0.0f);
        points[1] = EigenVector3<T>(1.0f, -1.0f, -1.0f)
                  + EigenVector3<T>(-0.5f, 0.0f, 0.0f);
        points[2] = EigenVector3<T>(1.0f, 1.0f, -1.0f)
                  + EigenVector3<T>(-0.5f, 0.0f, 0.0f);
        points[3] = EigenVector3<T>(-1.0f, 1.0f, -1.0f)
                  + EigenVector3<T>(-0.5f, 0.0f, 0.0f);
        points[4] = EigenVector3<T>(-1.0f, -1.0f, 1.0f)
                  + EigenVector3<T>(-0.5f, 0.0f, 0.0f);
        points[5] = EigenVector3<T>(1.0f, -1.0f, 1.0f)
                  + EigenVector3<T>(-0.5f, 0.0f, 0.0f);
        points[6] = EigenVector3<T>(1.0f, 1.0f, 1.0f)
                  + EigenVector3<T>(-0.5f, 0.0f, 0.0f);
        points[7] = EigenVector3<T>(-1.0f, 1.0f, 1.0f)
                  + EigenVector3<T>(-0.5f, 0.0f, 0.0f);

        geometry::DOP<T, 6> A = geometry::make_dop(&points[0], &points[0] + 8, DT);

        points[0] = EigenVector3<T>(-1.0f, -1.0f, -1.0f)
                  + EigenVector3<T>(0.5f, 0.0f, 0.0f);
        points[1] = EigenVector3<T>(1.0f, -1.0f, -1.0f)
                  + EigenVector3<T>(0.5f, 0.0f, 0.0f);
        points[2] = EigenVector3<T>(1.0f, 1.0f, -1.0f)
                  + EigenVector3<T>(0.5f, 0.0f, 0.0f);
        points[3] = EigenVector3<T>(-1.0f, 1.0f, -1.0f)
                  + EigenVector3<T>(0.5f, 0.0f, 0.0f);
        points[4] = EigenVector3<T>(-1.0f, -1.0f, 1.0f)
                  + EigenVector3<T>(0.5f, 0.0f, 0.0f);
        points[5] = EigenVector3<T>(1.0f, -1.0f, 1.0f)
                  + EigenVector3<T>(0.5f, 0.0f, 0.0f);
        points[6] = EigenVector3<T>(1.0f, 1.0f, 1.0f)
                  + EigenVector3<T>(0.5f, 0.0f, 0.0f);
        points[7] = EigenVector3<T>(-1.0f, 1.0f, 1.0f)
                  + EigenVector3<T>(0.5f, 0.0f, 0.0f);

        geometry::DOP<T, 6> B = geometry::make_dop(&points[0], &points[0] + 8, DT);

        BOOST_CHECK(geometry::overlap_dop_dop(A, B) == true);
        BOOST_CHECK(geometry::overlap_dop_dop(B, A) == true);
    }
  // Touching kDOPs
    {
        geometry::DirectionTable<T, 3> DT = geometry::make3<T>();

        std::vector<EigenVector3<T>> points(8);

        points[0] = EigenVector3<T>(-1.0f, -1.0f, -1.0f)
                  + EigenVector3<T>(-1.0f, 0.0f, 0.0f);
        points[1] = EigenVector3<T>(1.0f, -1.0f, -1.0f)
                  + EigenVector3<T>(-1.0f, 0.0f, 0.0f);
        points[2] = EigenVector3<T>(1.0f, 1.0f, -1.0f)
                  + EigenVector3<T>(-1.0f, 0.0f, 0.0f);
        points[3] = EigenVector3<T>(-1.0f, 1.0f, -1.0f)
                  + EigenVector3<T>(-1.0f, 0.0f, 0.0f);
        points[4] = EigenVector3<T>(-1.0f, -1.0f, 1.0f)
                  + EigenVector3<T>(-1.0f, 0.0f, 0.0f);
        points[5] = EigenVector3<T>(1.0f, -1.0f, 1.0f)
                  + EigenVector3<T>(-1.0f, 0.0f, 0.0f);
        points[6] = EigenVector3<T>(1.0f, 1.0f, 1.0f)
                  + EigenVector3<T>(-1.0f, 0.0f, 0.0f);
        points[7] = EigenVector3<T>(-1.0f, 1.0f, 1.0f)
                  + EigenVector3<T>(-1.0f, 0.0f, 0.0f);

        geometry::DOP<T, 6> A = geometry::make_dop(&points[0], &points[0] + 8, DT);

        points[0] = EigenVector3<T>(-1.0f, -1.0f, -1.0f)
                  + EigenVector3<T>(1.0f, 0.0f, 0.0f);
        points[1] = EigenVector3<T>(1.0f, -1.0f, -1.0f)
                  + EigenVector3<T>(1.0f, 0.0f, 0.0f);
        points[2] = EigenVector3<T>(1.0f, 1.0f, -1.0f)
                  + EigenVector3<T>(1.0f, 0.0f, 0.0f);
        points[3] = EigenVector3<T>(-1.0f, 1.0f, -1.0f)
                  + EigenVector3<T>(1.0f, 0.0f, 0.0f);
        points[4] = EigenVector3<T>(-1.0f, -1.0f, 1.0f)
                  + EigenVector3<T>(1.0f, 0.0f, 0.0f);
        points[5] = EigenVector3<T>(1.0f, -1.0f, 1.0f)
                  + EigenVector3<T>(1.0f, 0.0f, 0.0f);
        points[6] = EigenVector3<T>(1.0f, 1.0f, 1.0f)
                  + EigenVector3<T>(1.0f, 0.0f, 0.0f);
        points[7] = EigenVector3<T>(-1.0f, 1.0f, 1.0f)
                  + EigenVector3<T>(1.0f, 0.0f, 0.0f);

        geometry::DOP<T, 6> B = geometry::make_dop(&points[0], &points[0] + 8, DT);

        BOOST_CHECK(geometry::is_valid(A) == true);
        BOOST_CHECK(geometry::is_valid(B) == true);

        BOOST_CHECK(geometry::overlap_dop_dop(A, B) == true);
        BOOST_CHECK(geometry::overlap_dop_dop(B, A) == true);
    }
  // Identical kDOPs
    {
        geometry::DirectionTable<T, 3> DT = geometry::make3<T>();

        std::vector<EigenVector3<T>> points(8);

        points[0] = EigenVector3<T>(-1.0f, -1.0f, -1.0f);
        points[1] = EigenVector3<T>(1.0f, -1.0f, -1.0f);
        points[2] = EigenVector3<T>(1.0f, 1.0f, -1.0f);
        points[3] = EigenVector3<T>(-1.0f, 1.0f, -1.0f);
        points[4] = EigenVector3<T>(-1.0f, -1.0f, 1.0f);
        points[5] = EigenVector3<T>(1.0f, -1.0f, 1.0f);
        points[6] = EigenVector3<T>(1.0f, 1.0f, 1.0f);
        points[7] = EigenVector3<T>(-1.0f, 1.0f, 1.0f);

        geometry::DOP<T, 6> A = geometry::make_dop(&points[0], &points[0] + 8, DT);
        geometry::DOP<T, 6> B = geometry::make_dop(&points[0], &points[0] + 8, DT);

        BOOST_CHECK(geometry::overlap_dop_dop(A, B) == true);
        BOOST_CHECK(geometry::overlap_dop_dop(B, A) == true);
    }
}

BOOST_AUTO_TEST_SUITE_END();
