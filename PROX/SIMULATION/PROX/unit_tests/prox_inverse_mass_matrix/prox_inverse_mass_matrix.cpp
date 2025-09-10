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

    // Create a 6x6 sparse matrix (mass matrix block)
    Eigen::SparseMatrix<T> A(6, 6);

    // Fill the matrix with consecutive values starting from 0
    T value = 0;
    for (int i = 0; i < 6; ++i)
    {
        for (int j = 0; j < 6; ++j) { A.insert(i, j) = value++; }
    }

    // Check the initial values
    BOOST_CHECK_EQUAL(A.coeff(0, 0), T(0));
    BOOST_CHECK_EQUAL(A.coeff(0, 1), T(1));
    BOOST_CHECK_EQUAL(A.coeff(0, 2), T(2));
    BOOST_CHECK_EQUAL(A.coeff(0, 3), T(3));
    BOOST_CHECK_EQUAL(A.coeff(0, 4), T(4));
    BOOST_CHECK_EQUAL(A.coeff(0, 5), T(5));
    BOOST_CHECK_EQUAL(A.coeff(1, 0),
                      T(6)); // This is the 7th element in row-major order

    // Modify the (0,0) element
    A.coeffRef(0, 0) = T(2.5);

    // Create a copy
    Eigen::SparseMatrix<T> B = A;

    // Check the copied values
    BOOST_CHECK_EQUAL(B.coeff(0, 0), T(2.5));
    BOOST_CHECK_EQUAL(B.coeff(0, 1), T(1));
    BOOST_CHECK_EQUAL(B.coeff(0, 2), T(2));
    BOOST_CHECK_EQUAL(B.coeff(0, 3), T(3));
    BOOST_CHECK_EQUAL(B.coeff(0, 4), T(4));
    BOOST_CHECK_EQUAL(B.coeff(0, 5), T(5));
    BOOST_CHECK_EQUAL(B.coeff(1, 0), T(6));

    Eigen::Matrix<T, 6, 6> denseA = A.toDense();
    denseA = denseA.inverse();
    Eigen::SparseMatrix<T> A_inv = denseA.sparseView();

    // Check the inverted values (approximate due to floating point precision)
    const T tolerance = 1e-3;
    BOOST_CHECK_CLOSE(A_inv.coeff(0, 0), 1 / T(2.5), tolerance);
    BOOST_CHECK_CLOSE(A_inv.coeff(0, 1), 0, tolerance);
    BOOST_CHECK_CLOSE(A_inv.coeff(0, 2), -0.533333361, tolerance);
    BOOST_CHECK_CLOSE(A_inv.coeff(0, 3), 0.13333334, tolerance);
    BOOST_CHECK_CLOSE(A_inv.coeff(0, 4), -0.13333334, tolerance);
    BOOST_CHECK_CLOSE(A_inv.coeff(0, 5), 0.0666666701, tolerance);
    BOOST_CHECK_CLOSE(A_inv.coeff(1, 0), 0.352721691, tolerance);
}

BOOST_AUTO_TEST_SUITE_END();
