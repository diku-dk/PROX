#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/Dense>

#define BOOST_AUTO_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include <boost/test/unit_test_suite.hpp>
#include <boost/test/test_tools.hpp>

BOOST_AUTO_TEST_SUITE(prod_mass_jac);

using T = float;

// Helper function to fill a matrix with consecutive values starting from a given value
template <typename MatrixType>
void fill_matrix_consecutive(MatrixType& mat, T start_value)
{
    for (int i = 0; i < mat.rows(); ++i)
    {
        for (int j = 0; j < mat.cols(); ++j)
        {
            mat.insert(i, j) = start_value++;
        }
    }
}

// Helper function to fill a diagonal block matrix with consecutive values
void fill_diagonal_block(Eigen::SparseMatrix<T>& M, size_t block_idx,
                         T start_value)
{
    size_t start_row = block_idx * 6;
    size_t start_col = block_idx * 6;
    T value = start_value;

    for (int i = 0; i < 6; ++i)
    {
        for (int j = 0; j < 6; ++j)
        {
            M.insert(start_row + i, start_col + j) = value++;
        }
    }
}

BOOST_AUTO_TEST_CASE(mass_compress_prod_test_case)
{
    // Create mass matrix M (6x6)
    Eigen::SparseMatrix<T> M(6, 6);
    fill_diagonal_block(M, 0, 1.5f);
    M.coeffRef(0, 0) = 2.5f; // Override first element

    // Create transposed Jacobian JT (6x4)
    Eigen::SparseMatrix<T> JT(6, 4);
    fill_matrix_consecutive(JT, 2.0f);

    // Compute product C = M * JT
    Eigen::SparseMatrix<T> C = M * JT;

    // Check all values
    BOOST_CHECK_EQUAL(C.coeff(0, 0), 360);
    BOOST_CHECK_EQUAL(C.coeff(0, 1), 385);
    BOOST_CHECK_EQUAL(C.coeff(0, 2), 410);
    BOOST_CHECK_EQUAL(C.coeff(0, 3), 435);

    BOOST_CHECK_EQUAL(C.coeff(1, 0), 790);
    BOOST_CHECK_EQUAL(C.coeff(1, 1), 850);
    BOOST_CHECK_EQUAL(C.coeff(1, 2), 910);
    BOOST_CHECK_EQUAL(C.coeff(1, 3), 970);

    BOOST_CHECK_EQUAL(C.coeff(2, 0), 1222);
    BOOST_CHECK_EQUAL(C.coeff(2, 1), 1318);
    BOOST_CHECK_EQUAL(C.coeff(2, 2), 1414);
    BOOST_CHECK_EQUAL(C.coeff(2, 3), 1510);

    BOOST_CHECK_EQUAL(C.coeff(3, 0), 1654);
    BOOST_CHECK_EQUAL(C.coeff(3, 1), 1786);
    BOOST_CHECK_EQUAL(C.coeff(3, 2), 1918);
    BOOST_CHECK_EQUAL(C.coeff(3, 3), 2050);

    BOOST_CHECK_EQUAL(C.coeff(4, 0), 2086);
    BOOST_CHECK_EQUAL(C.coeff(4, 1), 2254);
    BOOST_CHECK_EQUAL(C.coeff(4, 2), 2422);
    BOOST_CHECK_EQUAL(C.coeff(4, 3), 2590);

    BOOST_CHECK_EQUAL(C.coeff(5, 0), 2518);
    BOOST_CHECK_EQUAL(C.coeff(5, 1), 2722);
    BOOST_CHECK_EQUAL(C.coeff(5, 2), 2926);
    BOOST_CHECK_EQUAL(C.coeff(5, 3), 3130);
}

BOOST_AUTO_TEST_CASE(mass_compress_resized_prod_test_case)
{
    // Create mass matrix M (12x12) with 2 blocks
    Eigen::SparseMatrix<T> M(12, 12);
    fill_diagonal_block(M, 0, 1.5f);
    M.coeffRef(0, 0) = 2.5f; // Override first element of first block

    fill_diagonal_block(M, 1, 2.5f); // Second block

    // Create transposed Jacobian JT (12x8) with 2x4 blocks
    Eigen::SparseMatrix<T> JT(12, 8);
    fill_matrix_consecutive(JT, 2.0f);

    // Set specific value for block (1,2)
    JT.coeffRef(6, 4) = 3.0f; // First element of block (1,2)

    // Compute product C = M * JT
    Eigen::SparseMatrix<T> C = M * JT;

    // Check values for block (0,0)
    BOOST_CHECK_EQUAL(C.coeff(0, 0), 670);
    BOOST_CHECK_EQUAL(C.coeff(0, 1), 695);
    BOOST_CHECK_EQUAL(C.coeff(0, 2), 720);
    BOOST_CHECK_EQUAL(C.coeff(0, 3), 745);
    BOOST_CHECK_EQUAL(C.coeff(1, 0), 1460);
    BOOST_CHECK_EQUAL(C.coeff(1, 1), 1520);
    BOOST_CHECK_EQUAL(C.coeff(1, 2), 1580);
    BOOST_CHECK_EQUAL(C.coeff(1, 3), 1640);
    BOOST_CHECK_EQUAL(C.coeff(2, 0), 2252);
    BOOST_CHECK_EQUAL(C.coeff(2, 1), 2348);
    BOOST_CHECK_EQUAL(C.coeff(2, 2), 2444);
    BOOST_CHECK_EQUAL(C.coeff(2, 3), 2540);
    BOOST_CHECK_EQUAL(C.coeff(3, 0), 3044);
    BOOST_CHECK_EQUAL(C.coeff(3, 1), 3176);
    BOOST_CHECK_EQUAL(C.coeff(3, 2), 3308);
    BOOST_CHECK_EQUAL(C.coeff(3, 3), 3440);
    BOOST_CHECK_EQUAL(C.coeff(4, 0), 3836);
    BOOST_CHECK_EQUAL(C.coeff(4, 1), 4004);
    BOOST_CHECK_EQUAL(C.coeff(4, 2), 4172);
    BOOST_CHECK_EQUAL(C.coeff(4, 3), 4340);
    BOOST_CHECK_EQUAL(C.coeff(5, 0), 4628);
    BOOST_CHECK_EQUAL(C.coeff(5, 1), 4832);
    BOOST_CHECK_EQUAL(C.coeff(5, 2), 5036);
    BOOST_CHECK_EQUAL(C.coeff(5, 3), 5240);

    // Check values for block (1,2)
    BOOST_CHECK_EQUAL(C.coeff(6, 4), 2232.5);
    BOOST_CHECK_EQUAL(C.coeff(6, 5), 2390);
    BOOST_CHECK_EQUAL(C.coeff(6, 6), 2420);
    BOOST_CHECK_EQUAL(C.coeff(6, 7), 2450);
    BOOST_CHECK_EQUAL(C.coeff(7, 4), 4590.5);
    BOOST_CHECK_EQUAL(C.coeff(7, 5), 5090);
    BOOST_CHECK_EQUAL(C.coeff(7, 6), 5156);
    BOOST_CHECK_EQUAL(C.coeff(7, 7), 5222);
    BOOST_CHECK_EQUAL(C.coeff(8, 4), 6948.5);
    BOOST_CHECK_EQUAL(C.coeff(8, 5), 7790);
    BOOST_CHECK_EQUAL(C.coeff(8, 6), 7892);
    BOOST_CHECK_EQUAL(C.coeff(8, 7), 7994);
    BOOST_CHECK_EQUAL(C.coeff(9, 4), 9306.5);
    BOOST_CHECK_EQUAL(C.coeff(9, 5), 10490);
    BOOST_CHECK_EQUAL(C.coeff(9, 6), 10628);
    BOOST_CHECK_EQUAL(C.coeff(9, 7), 10766);
    BOOST_CHECK_EQUAL(C.coeff(10, 4), 11664.5);
    BOOST_CHECK_EQUAL(C.coeff(10, 5), 13190);
    BOOST_CHECK_EQUAL(C.coeff(10, 6), 13364);
    BOOST_CHECK_EQUAL(C.coeff(10, 7), 13538);
    BOOST_CHECK_EQUAL(C.coeff(11, 4), 14022.5);
    BOOST_CHECK_EQUAL(C.coeff(11, 5), 15890);
    BOOST_CHECK_EQUAL(C.coeff(11, 6), 16100);
    BOOST_CHECK_EQUAL(C.coeff(11, 7), 16310);
}

BOOST_AUTO_TEST_CASE(compress_zero_init_clear_test_case)
{
    // Create a 6x4 sparse matrix
    Eigen::SparseMatrix<T> C(6, 4);

    // Check dimensions
    BOOST_CHECK_EQUAL(C.rows(), 6);
    BOOST_CHECK_EQUAL(C.cols(), 4);

    // Check that all elements are initially zero
    for (int i = 0; i < 6; ++i)
    {
        for (int j = 0; j < 4; ++j) { BOOST_CHECK_EQUAL(C.coeff(i, j), 0); }
    }

    // Clear the matrix by resizing to 0x0
    C.resize(0, 0);

    // Check that the matrix is now empty
    BOOST_CHECK_EQUAL(C.rows(), 0);
    BOOST_CHECK_EQUAL(C.cols(), 0);
    BOOST_CHECK_EQUAL(C.nonZeros(), 0);
}

BOOST_AUTO_TEST_SUITE_END();
