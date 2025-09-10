#include <sparse.h>
#include <prox_math_policy.h>

#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/Dense>

#define BOOST_AUTO_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include <boost/test/unit_test_suite.hpp>
#include <boost/test/test_tools.hpp>

BOOST_AUTO_TEST_SUITE(mass_block_jac);

BOOST_AUTO_TEST_CASE(mass_block_test_case)
{
    using T = float;

    // Create a 6x6 sparse matrix to represent the mass block
    Eigen::SparseMatrix<T> b(6, 6);

    // Initialize to zeros
    b.setZero();

    // Size checks
    BOOST_CHECK_EQUAL(b.rows(), 6);
    BOOST_CHECK_EQUAL(b.cols(), 6);
    BOOST_CHECK_EQUAL(b.nonZeros(), 0); // Initially no non-zero elements

    // Check if all elements are zero
    bool all_zeros = true;
    for (int i = 0; i < b.rows(); ++i)
    {
        for (int j = 0; j < b.cols(); ++j)
        {
            if (b.coeff(i, j) != 0)
            {
                all_zeros = false;
                break;
            }
        }
    }
    BOOST_CHECK(all_zeros);

    // Set element (0, 0) to 0.5
    b.coeffRef(0, 0) = 0.5f;
    BOOST_CHECK_EQUAL(b.coeff(0, 0), 0.5f);

    // Set element (1, 1) to 1.5
    b.coeffRef(1, 1) = 1.5f;
    BOOST_CHECK_EQUAL(b.coeff(1, 1), 1.5f);

    // Check that (2, 2) is still zero (not the same as (0, 0) or (1, 1))
    BOOST_CHECK_EQUAL(b.coeff(2, 2), 0.0f);

    // Create a copy
    Eigen::SparseMatrix<T> b_copy = b;

    // Check that all elements are equal
    for (int i = 0; i < b.rows(); ++i)
    {
        for (int j = 0; j < b.cols(); ++j)
        {
            BOOST_CHECK_EQUAL(b.coeff(i, j), b_copy.coeff(i, j));
        }
    }

    // Set additional elements in the copy
    b_copy.coeffRef(3, 3) = 0.5f;
    b_copy.coeffRef(4, 4) = 2.5f;
    b_copy.coeffRef(4, 5) = 3.5f;
    b_copy.coeffRef(3, 4) = 4.5f;
    b_copy.coeffRef(5, 4) = 5.5f;
    b_copy.coeffRef(3, 5) = 6.5f;
    b_copy.coeffRef(5, 5) = 7.5f;

    // Check symmetry
    BOOST_CHECK_EQUAL(b_copy.coeff(3, 4), 4.5f);
    BOOST_CHECK_EQUAL(b_copy.coeff(5, 4), 5.5f);
    BOOST_CHECK_EQUAL(b_copy.coeff(3, 5), 6.5f);

    // Assign back to b
    b = b_copy;

    // Check specific values
    BOOST_CHECK_EQUAL(b.coeff(1, 1), 1.5f);  // Was overwritten by the copy
    BOOST_CHECK_EQUAL(b.coeff(3, 3), 0.5f);  // Was overwritten by the copy
    BOOST_CHECK_EQUAL(b.coeff(3, 4), 4.5f);
    BOOST_CHECK_EQUAL(b.coeff(3, 5), 6.5f);
    BOOST_CHECK_EQUAL(b.coeff(4, 4), 2.5f);
    BOOST_CHECK_EQUAL(b.coeff(5, 4), 5.5f);
    BOOST_CHECK_EQUAL(b.coeff(5, 5), 7.5f);
}

BOOST_AUTO_TEST_SUITE_END();
