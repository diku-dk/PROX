#include "solvers/prox_gauss_seidel_solver.h"
#include "solvers/prox_jacobi_solver.h"
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include "eigenhelperfunctions.h"

#define BOOST_AUTO_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include <boost/test/unit_test_suite.hpp>
#include <boost/test/test_tools.hpp>

BOOST_AUTO_TEST_SUITE(math_policy_functions);

BOOST_AUTO_TEST_CASE(compute_z_test)
{
    using T = float;

    // Define Eigen types
    typedef Eigen::VectorX<T> vector_type;
    typedef std::vector<Eigen::Matrix<T, 4, 4>> diagonal4x4_type;
    typedef Eigen::SparseMatrix<T> compressed4x6_type;
    typedef Eigen::SparseMatrix<T> compressed6x4_type;

    // Initialize vectors and matrices
    vector_type x(4), b(4), z(4);
    diagonal4x4_type R(1);
    compressed4x6_type J(4, 6);
    compressed6x4_type WJT(6, 4);

    // Fill the vectors with specific values
    x << 13, 8, 12, 1;
    b << 16, 5, 9, 4;

    // Fill R (4x4 diagonal) with values starting from 1
    R[0].setZero();
    T value = 1.0f;
    for (int i = 0; i < 4; i++) { R[0](i, i) = value++; }

    // Fill J (4x6) with values starting from 3
    std::vector<Eigen::Triplet<T>> triplets;
    value = 3.0f;
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 6; j++)
        {
            triplets.push_back(Eigen::Triplet<T>(i, j, value++));
        }
    }
    J.setFromTriplets(triplets.begin(), triplets.end());
    triplets.clear();

    // Fill WJT (6x4) with values starting from 2
    value = 2.0f;
    for (int i = 0; i < 6; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            triplets.push_back(Eigen::Triplet<T>(i, j, value++));
        }
    }
    WJT.setFromTriplets(triplets.begin(), triplets.end());

    // Compute z = x - R(J*WJ^T *x + b) using the Eigen function
    prox::computeZ_Eigen(x, R, J, WJT, b, z);

    // Verify results
    BOOST_CHECK_CLOSE(z(0), -17002.0f, 1e-5f);
    BOOST_CHECK_CLOSE(z(1), -65896.0f, 1e-5f);
    BOOST_CHECK_CLOSE(z(2), -146700.0f, 1e-5f);
    BOOST_CHECK_CLOSE(z(3), -259387.0f, 1e-5f);

    // Change R and recalculate
    value = 2.0f;
    for (int i = 0; i < 4; i++) { R[0](i, i) = value++; }

    // Recompute z
    prox::computeZ_Eigen(x, R, J, WJT, b, z);

    // Verify new results
    BOOST_CHECK_CLOSE(z(0), -34017.0f, 1e-5f);
    BOOST_CHECK_CLOSE(z(1), -98848.0f, 1e-5f);
    BOOST_CHECK_CLOSE(z(2), -195604.0f, 1e-5f);
    BOOST_CHECK_CLOSE(z(3), -324234.0f, 1e-5f);

    // Change WJT and recalculate
    triplets.clear();
    value = 5.0f;
    for (int i = 0; i < 6; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            triplets.push_back(Eigen::Triplet<T>(i, j, value++));
        }
    }
    WJT.setFromTriplets(triplets.begin(), triplets.end());

    // Recompute z
    prox::computeZ_Eigen(x, R, J, WJT, b, z);

    BOOST_CHECK_CLOSE(z(0), -40749.0f, 1e-5f);
    BOOST_CHECK_CLOSE(z(1), -119962.0f, 1e-5f);
    BOOST_CHECK_CLOSE(z(2), -238444.0f, 1e-5f);
    BOOST_CHECK_CLOSE(z(3), -396144.0f, 1e-5f);
}

BOOST_AUTO_TEST_CASE(compute_norm_inf_test)
{
    using T = float;

    // Create an Eigen vector with initial size of 4
    Eigen::VectorX<T> x(4);

    // Fill the vector with values 0-3
    for (int i = 0; i < 4; i++) { x(i) = static_cast<T>(i); }

    // Compute the infinity norm
    T out = prox::computeInfNorm_Eigen(x);

    // Check that max value is 3
    BOOST_CHECK_CLOSE(out, 3.0f, 1e-5f);

    // Resize to check function works for multiple blocks (12 elements = 3 blocks of 4)
    x.resize(12);
    x.setZero(); // Initialize all elements to zero

    // Assign a large value to check that the multiple blocks are traversed
    // Set element at position 7 (which is 1*4 + 3 = 7) to 59
    x(7) = 59.0f;
    out = prox::computeInfNorm_Eigen(x);
    BOOST_CHECK_CLOSE(out, 59.0f, 1e-5f);

    // Assign a larger value to check that the multiple blocks are traversed
    // Set element at position 11 (which is 2*4 + 3 = 11) to 432
    x(11) = 432.0f;
    out = prox::computeInfNorm_Eigen(x);
    BOOST_CHECK_CLOSE(out, 432.0f, 1e-5f);

    // Assign a large negative value to check that absolute values are considered
    // Set element at position 3 (which is 0*4 + 3 = 3) to -1432
    x(3) = -1432.0f;
    out = prox::computeInfNorm_Eigen(x);
    BOOST_CHECK_CLOSE(out, 1432.0f, 1e-5f);
}

BOOST_AUTO_TEST_CASE(compute_wjt_test)
{
    // Define block types
    typedef Eigen::Matrix<float, 6, 6> Matrix6x6;
    typedef Eigen::Matrix<float, 4, 6> Matrix4x6;
    typedef Eigen::Matrix<float, 6, 4> Matrix6x4;

    // Create sparse matrices with appropriate dimensions
    // W: diagonal 6x6 block matrix (6x6 overall)
    Eigen::SparseMatrix<float> W(6, 6);

    // J: 4x6 block matrix (4x6 overall)
    Eigen::SparseMatrix<float> J(4, 6);

    // JT: 6x4 block matrix (6x4 overall)
    Eigen::SparseMatrix<float> JT(6, 4);

    // WJT: 6x4 block matrix (6x4 overall)
    Eigen::SparseMatrix<float> WJT(6, 4);

    // Fill the matrices with values using triplets
    std::vector<Eigen::Triplet<float>> triplets;

    // Fill W (6x6) with values starting from 1, incrementing by 1
    float value = 0.0f;
    for (int i = 0; i < 6; i++)
    {
        for (int j = 0; j < 6; j++)
        {
            triplets.push_back(Eigen::Triplet<float>(i, j, value++));
        }
    }
    W.setFromTriplets(triplets.begin(), triplets.end());
    triplets.clear();

    // Fill J (4x6) with values starting from 3, incrementing by 1
    value = 3.0f;
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 6; j++)
        {
            triplets.push_back(Eigen::Triplet<float>(i, j, value++));
        }
    }
    J.setFromTriplets(triplets.begin(), triplets.end());
    triplets.clear();

    // Fill JT (6x4) with values starting from 4, incrementing by 1
    value = 4.0f;
    for (int i = 0; i < 6; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            triplets.push_back(Eigen::Triplet<float>(i, j, value++));
        }
    }
    JT.setFromTriplets(triplets.begin(), triplets.end());

    // Compute WJT = W * J^T
    Eigen::SparseMatrix<float> J_transpose = J.transpose();
    WJT = W * J_transpose;

    // Convert to dense for easier verification
    Eigen::MatrixXf denseWJT = Eigen::MatrixXf(WJT);

    // Verify results - note the indexing difference
    BOOST_CHECK_CLOSE(denseWJT(0, 0), 100.0f, 1e-5f);
    BOOST_CHECK_CLOSE(denseWJT(0, 1), 190.0f, 1e-5f);
    BOOST_CHECK_CLOSE(denseWJT(0, 2), 280.0f, 1e-5f);
    BOOST_CHECK_CLOSE(denseWJT(0, 3), 370.0f, 1e-5f);
    BOOST_CHECK_CLOSE(denseWJT(1, 0), 298.0f, 1e-5f);
    BOOST_CHECK_CLOSE(denseWJT(1, 1), 604.0f, 1e-5f);
    BOOST_CHECK_CLOSE(denseWJT(1, 2), 910.0f, 1e-5f);
    BOOST_CHECK_CLOSE(denseWJT(1, 3), 1216.0f, 1e-5f);
    BOOST_CHECK_CLOSE(denseWJT(2, 0), 496.0f, 1e-5f);
    BOOST_CHECK_CLOSE(denseWJT(2, 1), 1018.0f, 1e-5f);
    BOOST_CHECK_CLOSE(denseWJT(2, 2), 1540.0f, 1e-5f);
    BOOST_CHECK_CLOSE(denseWJT(2, 3), 2062.0f, 1e-5f);
    BOOST_CHECK_CLOSE(denseWJT(3, 0), 694.0f, 1e-5f);
    BOOST_CHECK_CLOSE(denseWJT(3, 1), 1432.0f, 1e-5f);
    BOOST_CHECK_CLOSE(denseWJT(3, 2), 2170.0f, 1e-5f);
    BOOST_CHECK_CLOSE(denseWJT(3, 3), 2908.0f, 1e-5f);
    BOOST_CHECK_CLOSE(denseWJT(4, 0), 892.0f, 1e-5f);
    BOOST_CHECK_CLOSE(denseWJT(4, 1), 1846.0f, 1e-5f);
    BOOST_CHECK_CLOSE(denseWJT(4, 2), 2800.0f, 1e-5f);
    BOOST_CHECK_CLOSE(denseWJT(4, 3), 3754.0f, 1e-5f);
    BOOST_CHECK_CLOSE(denseWJT(5, 0), 1090.0f, 1e-5f);
    BOOST_CHECK_CLOSE(denseWJT(5, 1), 2260.0f, 1e-5f);
    BOOST_CHECK_CLOSE(denseWJT(5, 2), 3430.0f, 1e-5f);
    BOOST_CHECK_CLOSE(denseWJT(5, 3), 4600.0f, 1e-5f);

    // Second part: Compute WJT = W * JT
    WJT = W * JT;
    denseWJT = Eigen::MatrixXf(WJT);

    // Verify results
    BOOST_CHECK_CLOSE(denseWJT(0, 0), 280.0f, 1e-5f);
    BOOST_CHECK_CLOSE(denseWJT(0, 1), 295.0f, 1e-5f);
    BOOST_CHECK_CLOSE(denseWJT(0, 2), 310.0f, 1e-5f);
    BOOST_CHECK_CLOSE(denseWJT(0, 3), 325.0f, 1e-5f);
    BOOST_CHECK_CLOSE(denseWJT(1, 0), 784.0f, 1e-5f);
    BOOST_CHECK_CLOSE(denseWJT(1, 1), 835.0f, 1e-5f);
    BOOST_CHECK_CLOSE(denseWJT(1, 2), 886.0f, 1e-5f);
    BOOST_CHECK_CLOSE(denseWJT(1, 3), 937.0f, 1e-5f);
    BOOST_CHECK_CLOSE(denseWJT(2, 0), 1288.0f, 1e-5f);
    BOOST_CHECK_CLOSE(denseWJT(2, 1), 1375.0f, 1e-5f);
    BOOST_CHECK_CLOSE(denseWJT(2, 2), 1462.0f, 1e-5f);
    BOOST_CHECK_CLOSE(denseWJT(2, 3), 1549.0f, 1e-5f);
    BOOST_CHECK_CLOSE(denseWJT(3, 0), 1792.0f, 1e-5f);
    BOOST_CHECK_CLOSE(denseWJT(3, 1), 1915.0f, 1e-5f);
    BOOST_CHECK_CLOSE(denseWJT(3, 2), 2038.0f, 1e-5f);
    BOOST_CHECK_CLOSE(denseWJT(3, 3), 2161.0f, 1e-5f);
    BOOST_CHECK_CLOSE(denseWJT(4, 0), 2296.0f, 1e-5f);
    BOOST_CHECK_CLOSE(denseWJT(4, 1), 2455.0f, 1e-5f);
    BOOST_CHECK_CLOSE(denseWJT(4, 2), 2614.0f, 1e-5f);
    BOOST_CHECK_CLOSE(denseWJT(4, 3), 2773.0f, 1e-5f);
    BOOST_CHECK_CLOSE(denseWJT(5, 0), 2800.0f, 1e-5f);
    BOOST_CHECK_CLOSE(denseWJT(5, 1), 2995.0f, 1e-5f);
    BOOST_CHECK_CLOSE(denseWJT(5, 2), 3190.0f, 1e-5f);
    BOOST_CHECK_CLOSE(denseWJT(5, 3), 3385.0f, 1e-5f);
}

BOOST_AUTO_TEST_CASE(compute_b_test)
{
    using T = float;

    // Initialize vectors with appropriate sizes
    Eigen::VectorX<T> Wdth(6);
    Eigen::VectorX<T> u(6); // Same size as Wdth
    Eigen::VectorX<T> e(4); // 4 elements based on the checks
    Eigen::VectorX<T> g(4); // Same size as e
    Eigen::VectorX<T> b(4); // Output vector with 4 elements

    // Initialize sparse matrix J (4x6)
    Eigen::SparseMatrix<T> J(4, 6);

    // Fill J with incrementing values starting from 1.0f
    std::vector<Eigen::Triplet<T>> triplets;
    float value = 1.0f;
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 6; j++)
        {
            triplets.push_back(Eigen::Triplet<T>(i, j, value++));
        }
    }
    J.setFromTriplets(triplets.begin(), triplets.end());

    // Fill e with incrementing values starting from 2.0f
    value = 2.0f;
    for (int i = 0; i < 4; i++) { e(i) = value++; }

    // Verify e values
    BOOST_CHECK_CLOSE(e(0), 2.0f, 1e-5f);
    BOOST_CHECK_CLOSE(e(1), 3.0f, 1e-5f);
    BOOST_CHECK_CLOSE(e(2), 4.0f, 1e-5f);
    BOOST_CHECK_CLOSE(e(3), 5.0f, 1e-5f);

    // Fill Wdth with incrementing values starting from 1.0f
    value = 1.0f;
    for (int i = 0; i < 6; i++) { Wdth(i) = value++; }

    // Compute g = 2.0 * e and u = 2.0 * Wdth
    g = 2.0f * e;
    u = 2.0f * Wdth;

    // Call the Eigen version of compute_b
    prox::computeB_Eigen(J, Wdth, u, e, g, b);

    // Verify results (using CLOSE for floating point comparisons)
    BOOST_CHECK_CLOSE(b(0), 641.0f, 1e-5f);
    BOOST_CHECK_CLOSE(b(1), 1959.0f, 1e-5f);
    BOOST_CHECK_CLOSE(b(2), 3781.0f, 1e-5f);
    BOOST_CHECK_CLOSE(b(3), 6107.0f, 1e-5f);

    // Refill J with incrementing values starting from 3.0f
    triplets.clear();
    value = 3.0f;
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 6; j++)
        {
            triplets.push_back(Eigen::Triplet<T>(i, j, value++));
        }
    }
    J.setFromTriplets(triplets.begin(), triplets.end());

    // Call computeB_Eigen again with the updated J
    prox::computeB_Eigen(J, Wdth, u, e, g, b);

    // Verify the new results
    BOOST_CHECK_CLOSE(b(0), 935.0f, 1e-5f);
    BOOST_CHECK_CLOSE(b(1), 2337.0f, 1e-5f);
    BOOST_CHECK_CLOSE(b(2), 4243.0f, 1e-5f);
    BOOST_CHECK_CLOSE(b(3), 6653.0f, 1e-5f);
}

template <typename T>
Eigen::Matrix<T, 4, 1> get_block4(const Eigen::VectorX<T>& v, size_t k)
{
    return Eigen::Matrix<T, 4, 1>(v(4 * k), v(4 * k + 1), v(4 * k + 2),
                                  v(4 * k + 3));
}

template <typename T>
void computeZkA(const Eigen::Matrix<T, 4, 1>& x_k, const Eigen::VectorX<T>& w,
                const Eigen::Matrix<T, 4, 4>& R_k,
                const Eigen::SparseMatrix<T>& J,
                const Eigen::Matrix<T, 4, 1>& b_k, Eigen::Matrix<T, 4, 1>& z_k,
                const size_t& k)
{
    // Extract the k-th block row of J (rows 4*k to 4*k+3)
    Eigen::SparseMatrix<T> J_block = J;

    // Compute J_block * w
    z_k = J_block * w;

    // Add b_k
    z_k += b_k;

    // Multiply by R_k
    z_k = R_k * z_k;

    // Subtract from x_k
    z_k = x_k - z_k;
}

BOOST_AUTO_TEST_CASE(compute_z_k_test)
{
    // Note that the test is not a realistic usage example

    using T = float;
    Eigen::Matrix<T, 4, 4> R;
    Eigen::VectorX<T> Wdth(6);
    Eigen::VectorX<T> w(6); // Same size as Wdth
    Eigen::VectorX<T> x(4);
    Eigen::VectorX<T> b(4);
    Eigen::Matrix<T, 4, 1> z_k;

    // Initialize sparse matrix J (4x6)
    Eigen::SparseMatrix<T> J(4, 6);
    Eigen::SparseMatrix<T> WJT(6, 4);

    // Fill the vectors with specific values
    x << 13, 8, 12, 1;
    b << 16, 5, 9, 4;

    // Fill the matrices with incrementing values
    std::vector<Eigen::Triplet<float>> triplets;

    R.setZero();
    T value = 1.0f;
    for (int i = 0; i < 4; i++) { R(i, i) = value++; }

    // Fill J (4x6) with values starting from 3
    value = 3.0f;
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 6; j++)
        {
            triplets.push_back(Eigen::Triplet<float>(i, j, value++));
        }
    }
    J.setFromTriplets(triplets.begin(), triplets.end());
    triplets.clear();

    // Fill WJT (6x4) with values starting from 2
    value = 2.0f;
    for (int i = 0; i < 6; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            triplets.push_back(Eigen::Triplet<float>(i, j, value++));
        }
    }
    WJT.setFromTriplets(triplets.begin(), triplets.end());

    // Compute w = WJT * x
    w = WJT * x;

    // Compute z_k using the Eigen function
    prox::computeZk_Eigen(get_block4(x, 0), w, R, J, get_block4(b, 0), z_k, 0);

    // Verify results
    BOOST_CHECK_CLOSE(z_k(0), -17002.0f, 1e-5f);
    BOOST_CHECK_CLOSE(z_k(1), -65896.0f, 1e-5f);
    BOOST_CHECK_CLOSE(z_k(2), -146700.0f, 1e-5f);
    BOOST_CHECK_CLOSE(z_k(3), -259387.0f, 1e-5f);

    // Change R and recalculate
    value = 2.0f;
    for (int i = 0; i < 4; i++) { R(i, i) = value++; }

    // Recompute w = WJT * x
    w = WJT * x;

    // Recompute z_k
    prox::computeZk_Eigen(get_block4(x, 0), w, R, J, get_block4(b, 0), z_k, 0);

    // Verify new results
    BOOST_CHECK_CLOSE(z_k(0), -34017.0f, 1e-5f);
    BOOST_CHECK_CLOSE(z_k(1), -98848.0f, 1e-5f);
    BOOST_CHECK_CLOSE(z_k(2), -195604.0f, 1e-5f);
    BOOST_CHECK_CLOSE(z_k(3), -324234.0f, 1e-5f);

    // Change WJT and recalculate
    triplets.clear();
    value = 5.0f;
    for (int i = 0; i < 6; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            triplets.push_back(Eigen::Triplet<float>(i, j, value++));
        }
    }
    WJT.setFromTriplets(triplets.begin(), triplets.end());

    // Recompute w = WJT * x
    w = WJT * x;

    // Recompute z_k
    prox::computeZk_Eigen(get_block4(x, 0), w, R, J, get_block4(b, 0), z_k, 0);

    // Verify new results
    BOOST_CHECK_CLOSE(z_k(0), -40749.0f, 1e-5f);
    BOOST_CHECK_CLOSE(z_k(1), -119962.0f, 1e-5f);
    BOOST_CHECK_CLOSE(z_k(2), -238444.0f, 1e-5f);
    BOOST_CHECK_CLOSE(z_k(3), -396144.0f, 1e-5f);
}

BOOST_AUTO_TEST_SUITE_END();
