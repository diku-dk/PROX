#include <geometry.h>
#include <eigenhelperall.h>

#define BOOST_AUTO_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include <boost/test/unit_test_suite.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/test/test_tools.hpp>

using T = float;

BOOST_AUTO_TEST_SUITE(geometry);

BOOST_AUTO_TEST_CASE(overlap_tri_tri_test)
{
  // B inside A
    {
        EigenVector3<T> const A0 = EigenVector3<T>(0.0, 0.0, 0.0);
        EigenVector3<T> const A1 = EigenVector3<T>(1.0, 0.0, 0.0);
        EigenVector3<T> const A2 = EigenVector3<T>(0.0, 1.0, 0.0);

        EigenVector3<T> const B0 = EigenVector3<T>(0.1, 0.1, 0.1);
        EigenVector3<T> const B1 = EigenVector3<T>(0.9, 0.0, 0.0);
        EigenVector3<T> const B2 = EigenVector3<T>(0.0, 0.9, 0.0);

        geometry::Triangle<T> triA = geometry::make_triangle<T>(A0, A1, A2);
        geometry::Triangle<T> triB = geometry::make_triangle<T>(B0, B1, B2);

        bool const test1 = geometry::overlap_triangle_triangle(triA, triB);
        bool const test2 = geometry::overlap_triangle_triangle(triB, triA);

        BOOST_CHECK(test1);
        BOOST_CHECK(test2);
    }

    unsigned int permutation[6][3] = {
        {0, 1, 2},
        {0, 2, 1},
        {1, 2, 0},
        {1, 0, 2},
        {2, 0, 1},
        {2, 1, 0}
    };

  // Touching edge-edge cases
    {
        std::vector<EigenVector3<T>> A(3u);
        std::vector<EigenVector3<T>> B(3u);

        A[0] = EigenVector3<T>(0.0, 0.0, 0.0);
        A[1] = EigenVector3<T>(1.0, 0.0, 0.0);
        A[2] = EigenVector3<T>(0.0, 1.0, 0.0);

        B[0] = EigenVector3<T>(0.5, 0.0, -0.5);
        B[1] = EigenVector3<T>(0.5, 0.0, 0.5);
        B[2] = EigenVector3<T>(0.5, -1.0, 0.0);

        for (unsigned int i = 0u; i < 6u; ++i)
        {
            for (unsigned int j = 0u; j < 6u; ++j)
            {
                geometry::Triangle<T> triA = geometry::make_triangle<T>(
                    A[permutation[i][0]], A[permutation[i][1]],
                    A[permutation[i][2]]);
                geometry::Triangle<T> triB = geometry::make_triangle<T>(
                    B[permutation[j][0]], B[permutation[j][1]],
                    B[permutation[j][2]]);

                bool const test1 = geometry::overlap_triangle_triangle(triA, triB);
                bool const test2 = geometry::overlap_triangle_triangle(triB, triA);

                BOOST_CHECK(test1);
                BOOST_CHECK(test2);
            }
        }
    }

  // Separating edge-edge cases
    {
        std::vector<EigenVector3<T>> A(3u);
        std::vector<EigenVector3<T>> B(3u);

        A[0] = EigenVector3<T>(0.0, 0.0, 0.0);
        A[1] = EigenVector3<T>(1.0, 0.0, 0.0);
        A[2] = EigenVector3<T>(0.0, 1.0, 0.0);

        B[0] = EigenVector3<T>(0.5, -0.01, -0.5);
        B[1] = EigenVector3<T>(0.5, -0.01, 0.5);
        B[2] = EigenVector3<T>(0.5, -1.01, 0.0);

        for (unsigned int i = 0u; i < 6u; ++i)
        {
            for (unsigned int j = 0u; j < 6u; ++j)
            {
                geometry::Triangle<T> triA = geometry::make_triangle<T>(
                    A[permutation[i][0]], A[permutation[i][1]],
                    A[permutation[i][2]]);
                geometry::Triangle<T> triB = geometry::make_triangle<T>(
                    B[permutation[j][0]], B[permutation[j][1]],
                    B[permutation[j][2]]);

                bool const test1 = geometry::overlap_triangle_triangle(triA, triB);
                bool const test2 = geometry::overlap_triangle_triangle(triB, triA);

                BOOST_CHECK(!test1);
                BOOST_CHECK(!test2);
            }
        }
    }

  // Separated by face cases
    {
        std::vector<EigenVector3<T>> A(3u);
        std::vector<EigenVector3<T>> B(3u);

        A[0] = EigenVector3<T>(0.0, 0.0, 0.0);
        A[1] = EigenVector3<T>(1.0, 0.0, 0.0);
        A[2] = EigenVector3<T>(0.0, 1.0, 0.0);

        B[0] = EigenVector3<T>(1.1, 0.0, 0.1);
        B[1] = EigenVector3<T>(1.1, 0.0, 1.1);
        B[2] = EigenVector3<T>(1.1, -1.0, 0.6);

        for (unsigned int i = 0u; i < 6u; ++i)
        {
            for (unsigned int j = 0u; j < 6u; ++j)
            {
                geometry::Triangle<T> triA = geometry::make_triangle<T>(
                    A[permutation[i][0]], A[permutation[i][1]],
                    A[permutation[i][2]]);
                geometry::Triangle<T> triB = geometry::make_triangle<T>(
                    B[permutation[j][0]], B[permutation[j][1]],
                    B[permutation[j][2]]);

                bool const test1 = geometry::overlap_triangle_triangle(triA, triB);
                bool const test2 = geometry::overlap_triangle_triangle(triB, triA);

                BOOST_CHECK(!test1);
                BOOST_CHECK(!test2);
            }
        }
    }
}

BOOST_AUTO_TEST_SUITE_END();
