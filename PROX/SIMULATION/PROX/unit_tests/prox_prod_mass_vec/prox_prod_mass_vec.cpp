#include <sparse.h>
#include <sparse_fill.h>
#include <prox_math_policy.h>

#define BOOST_AUTO_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include <boost/test/unit_test_suite.hpp>
#include <boost/test/test_tools.hpp>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>

BOOST_AUTO_TEST_SUITE(prod_mass_vec);

// Helper function to fill a diagonal block matrix
template <typename T>
void fill_diagonal_block(Eigen::SparseMatrix<T>& M, size_t block_idx,
                         const std::vector<T>& values)
{
    size_t start = block_idx * 6;
    for (int i = 0; i < 6; ++i) { M.insert(start + i, start + i) = values[i]; }
}

// Helper function to fill a vector block
template <typename T>
void fill_vector_block(Eigen::VectorX<T>& v, size_t block_idx, T start_value)
{
    size_t start = block_idx * 6;
    for (int i = 0; i < 6; ++i) { v(start + i) = start_value + i; }
}

BOOST_AUTO_TEST_CASE(mass_vec_prod_test_case)
{
    using T = float;

    // Test 1: Single block
    {
        Eigen::SparseMatrix<T> M(6, 6);
        std::vector<T> diag_values = {2.5, 1, 2, 3, 4, 5};
        fill_diagonal_block(M, 0, diag_values);

        Eigen::VectorX<T> u(6);
        fill_vector_block(u, 0, 1.0f);

        Eigen::VectorX<T> r = M * u;

        BOOST_CHECK_EQUAL(r(0), T(2.5));
        BOOST_CHECK_EQUAL(r(1), T(2));
        BOOST_CHECK_EQUAL(r(2), T(6));
        BOOST_CHECK_EQUAL(r(3), T(12));
        BOOST_CHECK_EQUAL(r(4), T(20));
        BOOST_CHECK_EQUAL(r(5), T(30));
    }

    // Test 2: Two blocks
    {
        Eigen::SparseMatrix<T> M(12, 12);

        // First block
        std::vector<T> diag_values1 = {2.5, 1, 2, 3, 4, 5};
        fill_diagonal_block(M, 0, diag_values1);

        // Second block
        std::vector<T> diag_values2 = {2.5, 1, 2, 3, 4, 5};
        fill_diagonal_block(M, 1, diag_values2);

        Eigen::VectorX<T> u(12);
        fill_vector_block(u, 0, 1.0f);
        fill_vector_block(u, 1, 7.0f);

        Eigen::VectorX<T> r = M * u;

        // First block results
        BOOST_CHECK_EQUAL(r(0), T(2.5));
        BOOST_CHECK_EQUAL(r(1), T(2));
        BOOST_CHECK_EQUAL(r(2), T(6));
        BOOST_CHECK_EQUAL(r(3), T(12));
        BOOST_CHECK_EQUAL(r(4), T(20));
        BOOST_CHECK_EQUAL(r(5), T(30));

        // Second block results
        BOOST_CHECK_EQUAL(r(6), T(17.5));
        BOOST_CHECK_EQUAL(r(7), T(8));
        BOOST_CHECK_EQUAL(r(8), T(18));
        BOOST_CHECK_EQUAL(r(9), T(30));
        BOOST_CHECK_EQUAL(r(10), T(44));
        BOOST_CHECK_EQUAL(r(11), T(60));
    }
}

BOOST_AUTO_TEST_SUITE_END();
