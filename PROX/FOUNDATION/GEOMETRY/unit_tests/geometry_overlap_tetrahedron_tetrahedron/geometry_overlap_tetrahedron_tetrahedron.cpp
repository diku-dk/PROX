#include <geometry.h>
#include <tiny.h>

#define BOOST_AUTO_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include <boost/test/unit_test_suite.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/test/test_tools.hpp>

using T = float;

BOOST_AUTO_TEST_SUITE(geometry);

BOOST_AUTO_TEST_CASE(overlap_tet_tet_test)
{
  // B inside A, all SAT tests should fail
    {
        EigenVector3<T> const Ap0 = EigenVector3<T>(0.0, 0.0, 0.0);
        EigenVector3<T> const Ap1 = EigenVector3<T>(1.0, 0.0, 0.0);
        EigenVector3<T> const Ap2 = EigenVector3<T>(0.0, 1.0, 0.0);
        EigenVector3<T> const Ap3 = EigenVector3<T>(0.0, 0.0, 1.0);

        EigenVector3<T> const Bp0 = EigenVector3<T>(0.1, 0.1, 0.1);
        EigenVector3<T> const Bp1 = EigenVector3<T>(0.9, 0.0, 0.0);
        EigenVector3<T> const Bp2 = EigenVector3<T>(0.0, 0.9, 0.0);
        EigenVector3<T> const Bp3 = EigenVector3<T>(0.0, 0.0, 0.9);

        geometry::TetrahedronEigen<T> tetA
            = geometry::make_tetrahedron<T>(Ap0, Ap1, Ap2, Ap3);
        geometry::TetrahedronEigen<T> tetB
            = geometry::make_tetrahedron<T>(Bp0, Bp1, Bp2, Bp3);

        bool const test1 = geometry::overlap_tetrahedron_tetrahedron(tetA, tetB);
        bool const test2 = geometry::overlap_tetrahedron_tetrahedron(tetB, tetA);

        BOOST_CHECK(test1);
        BOOST_CHECK(test2);
    }

    unsigned int permutation[6][4] = {
        {1, 3, 2, 0},
        {0, 3, 1, 2},
        {0, 1, 2, 3},
        {2, 3, 0, 1},
        {1, 2, 0, 3},
        {2, 0, 1, 3}
    };

  // Touching edge-edge cases
    {
        std::vector<EigenVector3<T>> A(4u);
        std::vector<EigenVector3<T>> B(4u);

        A[0] = EigenVector3<T>(0.0, 0.0, 0.0);
        A[1] = EigenVector3<T>(1.0, 0.0, 0.0);
        A[2] = EigenVector3<T>(0.0, 1.0, 0.0);
        A[3] = EigenVector3<T>(0.0, 0.0, 1.0);

        B[0] = EigenVector3<T>(1.0, -1.0, 0.5);
        B[1] = EigenVector3<T>(-1.0, 1.0, 0.5);
        B[2] = EigenVector3<T>(-1.0, -1.0, -1.0);
        B[3] = EigenVector3<T>(-1.0, -1.0, 1.0);

        for (unsigned int i = 0u; i < 6u; ++i)
        {
            for (unsigned int j = 0u; j < 6u; ++j)
            {

                geometry::TetrahedronEigen<T> tetA
                    = geometry::make_tetrahedron<T>(
                        A[permutation[i][0]], A[permutation[i][1]],
                        A[permutation[i][2]], A[permutation[i][3]]);
                geometry::TetrahedronEigen<T> tetB
                    = geometry::make_tetrahedron<T>(
                        B[permutation[j][0]], B[permutation[j][1]],
                        B[permutation[j][2]], B[permutation[j][3]]);

                bool const test1 = geometry::overlap_tetrahedron_tetrahedron(tetA, tetB);
                bool const test2 = geometry::overlap_tetrahedron_tetrahedron(tetB, tetA);

                BOOST_CHECK(test1);
                BOOST_CHECK(test2);
            }
        }
    }

  // Separating edge-edge cases
    {
        std::vector<EigenVector3<T>> A(4u);
        std::vector<EigenVector3<T>> B(4u);

        A[0] = EigenVector3<T>(0.01, 0.01, 0.01);
        A[1] = EigenVector3<T>(1.01, 0.01, 0.01);
        A[2] = EigenVector3<T>(0.01, 1.01, 0.01);
        A[3] = EigenVector3<T>(0.01, 0.01, 1.01);

        B[0] = EigenVector3<T>(1.0, -1.0, 0.5);
        B[1] = EigenVector3<T>(-1.0, 1.0, 0.5);
        B[2] = EigenVector3<T>(-1.0, -1.0, -1.0);
        B[3] = EigenVector3<T>(-1.0, -1.0, 1.0);

        for (unsigned int i = 0u; i < 6u; ++i)
        {
            for (unsigned int j = 0u; j < 6u; ++j)
            {

                geometry::TetrahedronEigen<T> tetA
                    = geometry::make_tetrahedron<T>(
                        A[permutation[i][0]], A[permutation[i][1]],
                        A[permutation[i][2]], A[permutation[i][3]]);
                geometry::TetrahedronEigen<T> tetB
                    = geometry::make_tetrahedron<T>(
                        B[permutation[j][0]], B[permutation[j][1]],
                        B[permutation[j][2]], B[permutation[j][3]]);

                bool const test1 = geometry::overlap_tetrahedron_tetrahedron(tetA, tetB);
                bool const test2 = geometry::overlap_tetrahedron_tetrahedron(tetB, tetA);

                BOOST_CHECK(!test1);
                BOOST_CHECK(!test2);
            }
        }
    }

  // Separated by face cases
    {
        EigenVector3<T> const A0 = EigenVector3<T>(0.0, 0.0, 0.0);
        EigenVector3<T> const A1 = EigenVector3<T>(1.0, 0.0, 0.0);
        EigenVector3<T> const A2 = EigenVector3<T>(0.0, 1.0, 0.0);
        EigenVector3<T> const A3 = EigenVector3<T>(0.0, 0.0, 1.0);

        std::vector<EigenVector3<T>> offset(4);

        offset[0] = EigenVector3<T>(2.0, 2.0, 2.0);
        offset[1] = EigenVector3<T>(-2.0, 0.0, 0.0);
        offset[2] = EigenVector3<T>(0.0, -2.0, 0.0);
        offset[3] = EigenVector3<T>(0.0, 0.0, -2.0);

        geometry::TetrahedronEigen<T> tetA
            = geometry::make_tetrahedron<T>(A0, A1, A2, A3);

        for (unsigned int k = 0u; k < 4u; ++k)
        {

            geometry::TetrahedronEigen<T> tetB = geometry::make_tetrahedron<T>(
                A0 + offset[k], A1 + offset[k], A2 + offset[k], A3 + offset[k]);

            bool const test1 = geometry::overlap_tetrahedron_tetrahedron(tetA, tetB);
            bool const test2 = geometry::overlap_tetrahedron_tetrahedron(tetB, tetA);

            BOOST_CHECK(!test1);
            BOOST_CHECK(!test2);
        }
    }
}

BOOST_AUTO_TEST_SUITE_END();
