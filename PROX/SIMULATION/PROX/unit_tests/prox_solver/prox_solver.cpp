#include <solvers/prox_jacobi_solver.h>
#include <solvers/prox_gauss_seidel_solver.h>

#define BOOST_AUTO_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include <boost/test/unit_test_suite.hpp>
#include <boost/test/test_tools.hpp>

/**
 * This file contain a number of tests for sequential prox solvers.
 *
 * The first three functions take a pointer to a prox solver function, setup a test
 * case and call the specified solver. The test functions can be used to unit test
 * a new solver by making a math_policy and calling the test functions.
 *
 * Following the test setup functions, the solvers are tested. As Boost only specify
 * the line and function where an error occurs, it is suggested that a test function is
 * made for each (solver, setup function) pair.
 */

/**
 * This function build a simple test setup for use with a solver for the contact problem,
 * and performs a number of tests. It takes the math policy and the solver function as
 * template arguments.
 *
 * The setup has a unit sphere placed at (0,2,0), resting
 * on top of a 2x2x2 fixed box with its center at the origin.
 * The expected result is that only the normal impluse is different from zero.
 * The speed of the sphere is g*t = -9.81*0.01 = -0.0981, and the impulse is
 * P = m*v = 20*-0.0981 = -1.962, so the magnitude of the normal impulse must
 * be close to 1.962.
 */
/*template <typename math_policy, typename prox_solver_functor>
void run_2_body_prox_solver_test(prox_solver_functor const& prox_solver)
{

    typename math_policy::vector7_type q(2);
    typename math_policy::vector6_type Wdth(2), u(2);
    typename math_policy::vector4_type lambda(1), g(1), e(1), mu(1);
    typename math_policy::compressed4x6_type J(1, 2, 2);
    typename math_policy::compressed6x4_type WJT(2, 1, 2);

    prox::SolverParams<math_policy> params;

  //filling the vectors
    q(0)[0] = 0;
    q(0)[1] = 0;
    q(0)[2] = 0;
    q(0)[3] = 0;
    q(0)[4] = 0;
    q(0)[5] = 0;
    q(0)[6] = -1.0;
    q(1)[0] = 0;
    q(1)[1] = 2.0;
    q(1)[2] = 0;
    q(1)[3] = 0;
    q(1)[4] = 0;
    q(1)[5] = 0;
    q(1)[6] = -1.0;

    u(0)[0] = 0;
    u(0)[1] = 0;
    u(0)[2] = 0;
    u(0)[3] = 0;
    u(0)[4] = 0;
    u(0)[5] = 0;
    u(1)[0] = 0;
    u(1)[1] = 0;
    u(1)[2] = 0;
    u(1)[3] = 0;
    u(1)[4] = 0;
    u(1)[5] = 0;

    Wdth(0)[0] = 0;
    Wdth(0)[1] = 0;
    Wdth(0)[2] = 0;
    Wdth(0)[3] = 0;
    Wdth(0)[4] = 0;
    Wdth(0)[5] = 0;
    Wdth(1)[0] = 0;
    Wdth(1)[1] = -0.0981;
    Wdth(1)[2] = 0;
    Wdth(1)[3] = 0;
    Wdth(1)[4] = 0;
    Wdth(1)[5] = 0;

    g(0)[0] = 0;
    g(0)[1] = 0;
    g(0)[2] = 0;
    g(0)[3] = 0;
    e(0)[0] = 1;
    e(0)[1] = 0;
    e(0)[2] = 0;
    e(0)[3] = 0;
    mu(0)[0] = 0;
    mu(0)[1] = 0.5;
    mu(0)[2] = 0.5;
    mu(0)[3] = 0.5;

  //filling the matrices
  //math_policy::block4x6_type & Jblock = J(0,0);
    J(0, 0)[0] = 0;
    J(0, 0)[1] = -1;
    J(0, 0)[2] = 0;
    J(0, 0)[3] = 0;
    J(0, 0)[4] = 0;
    J(0, 0)[5] = 0;
    J(0, 0)[6] = -1;
    J(0, 0)[7] = 0;
    J(0, 0)[8] = 0;
    J(0, 0)[9] = 0;
    J(0, 0)[10] = 0;
    J(0, 0)[11] = 1;
    J(0, 0)[12] = 0;
    J(0, 0)[13] = 0;
    J(0, 0)[14] = -1;
    J(0, 0)[15] = -1;
    J(0, 0)[16] = 0;
    J(0, 0)[17] = 0;
    J(0, 0)[18] = 0;
    J(0, 0)[19] = 0;
    J(0, 0)[20] = 0;
    J(0, 0)[21] = 0;
    J(0, 0)[22] = -1;
    J(0, 0)[23] = 0;

  //math_policy::block4x6_type & Jblock2 = J(0,1);
    J(0, 1)[0] = 0;
    J(0, 1)[1] = 1;
    J(0, 1)[2] = 0;
    J(0, 1)[3] = 0;
    J(0, 1)[4] = 0;
    J(0, 1)[5] = 0;
    J(0, 1)[6] = 1;
    J(0, 1)[7] = 0;
    J(0, 1)[8] = 0;
    J(0, 1)[9] = 0;
    J(0, 1)[10] = 0;
    J(0, 1)[11] = 1;
    J(0, 1)[12] = 0;
    J(0, 1)[13] = 0;
    J(0, 1)[14] = 1;
    J(0, 1)[15] = -1;
    J(0, 1)[16] = 0;
    J(0, 1)[17] = 0;
    J(0, 1)[18] = 0;
    J(0, 1)[19] = 0;
    J(0, 1)[20] = 0;
    J(0, 1)[21] = 0;
    J(0, 1)[22] = 1;
    J(0, 1)[23] = 0;

  //math_policy::block6x4_type & WJTblock = WJT(0,0);
    WJT(0, 0)[0] = 0;
    WJT(0, 0)[1] = 0;
    WJT(0, 0)[2] = 0;
    WJT(0, 0)[3] = 0;
    WJT(0, 0)[4] = 0;
    WJT(0, 0)[5] = 0;
    WJT(0, 0)[6] = 0;
    WJT(0, 0)[7] = 0;
    WJT(0, 0)[8] = 0;
    WJT(0, 0)[9] = 0;
    WJT(0, 0)[10] = 0;
    WJT(0, 0)[11] = 0;
    WJT(0, 0)[12] = 0;
    WJT(0, 0)[13] = 0;
    WJT(0, 0)[14] = 0;
    WJT(0, 0)[15] = 0;
    WJT(0, 0)[16] = 0;
    WJT(0, 0)[17] = 0;
    WJT(0, 0)[18] = 0;
    WJT(0, 0)[19] = 0;
    WJT(0, 0)[20] = 0;
    WJT(0, 0)[21] = 0;
    WJT(0, 0)[22] = 0;
    WJT(0, 0)[23] = 0;
  //WJTblock = WJT(1,0);
    WJT(1, 0)[0] = 0;
    WJT(1, 0)[1] = 0.05;
    WJT(1, 0)[2] = 0;
    WJT(1, 0)[3] = 0;
    WJT(1, 0)[4] = 0.05;
    WJT(1, 0)[5] = 0;
    WJT(1, 0)[6] = 0;
    WJT(1, 0)[7] = 0;
    WJT(1, 0)[8] = 0;
    WJT(1, 0)[9] = 0;
    WJT(1, 0)[10] = 0.05;
    WJT(1, 0)[11] = 0;
    WJT(1, 0)[12] = 0;
    WJT(1, 0)[13] = 0;
    WJT(1, 0)[14] = -0.05;
    WJT(1, 0)[15] = 0;
    WJT(1, 0)[16] = 0;
    WJT(1, 0)[17] = 0;
    WJT(1, 0)[18] = 0;
    WJT(1, 0)[19] = 0.05;
    WJT(1, 0)[20] = 0;
    WJT(1, 0)[21] = 0.05;
    WJT(1, 0)[22] = 0;
    WJT(1, 0)[23] = 0;

    BOOST_CHECK(false);
  //  prox_solver( ... fix me .... );
    typename math_policy::real_type tolerance = 1e-1; // that is 0.1 percent, an error on the fourth decimal
    BOOST_CHECK_CLOSE(lambda(0)[0], 1.962f, tolerance);
    BOOST_CHECK_EQUAL(lambda(0)[1], 0);
    BOOST_CHECK_EQUAL(lambda(0)[2], 0);
    BOOST_CHECK_EQUAL(lambda(0)[3], 0);

  // We now apply a small velocity in the direction of the x axis on the sphere.
  // This should result in friction. If the speed is 0.25, the impulse of the sphere
  // is P = m*v = 20*0.25 = 5. The bound on the friction impulse is
  // mu*l_n = 0.5*1.962 = 0.981.
    lambda.clear();
    u(0)[0] = 0.25;
    BOOST_CHECK(false);
  //  prox_solver( ... fix me .... );

    BOOST_CHECK_CLOSE(lambda(0)[0], 1.962f, tolerance);
    BOOST_CHECK_CLOSE(lambda(0)[1], 0.981f, tolerance);
    BOOST_CHECK_EQUAL(lambda(0)[2], 0);
    BOOST_CHECK_EQUAL(lambda(0)[3], 0);

  // Applying a small velocity in the direction of the z axis on the sphere, should
  // result in friction. If the speed is 0.01, the impulse is P = m*v = 20*0.01 = 0.2.
  // This changes the bounds on the friction to (0.980, 0.0390)
    lambda.clear();
    u(0)[2] = 0.01;

    BOOST_CHECK(false);
  //  prox_solver( ... fix me .... );

    BOOST_CHECK_CLOSE(lambda(0)[0], 1.962f, tolerance);
    BOOST_CHECK_CLOSE(lambda(0)[1], 0.980f, tolerance);
    BOOST_CHECK_CLOSE(lambda(0)[2], 0.0392f, tolerance);
    BOOST_CHECK_EQUAL(lambda(0)[3], 0);

  // We now apply an angular velocity around the y axis on the sphere, to test Contensou friction
  // Setting w_y = 0.25, yields an upper bound of mu l_n = 0.981 (we assume the 2/3 factor to be
  // included in mu)
    u(0)[0] = 0;
    u(0)[2] = 0;
    u(0)[4] = 0.25;
    lambda.clear();

    BOOST_CHECK(false);
  //  prox_solver( ... fix me .... );

    BOOST_CHECK_CLOSE(lambda(0)[0], 1.962f, tolerance);
    BOOST_CHECK_EQUAL(lambda(0)[1], 0);
    BOOST_CHECK_EQUAL(lambda(0)[2], 0);
    BOOST_CHECK_CLOSE(lambda(0)[3], 0.981f, tolerance);
}*/
template <typename T>
void set_block6(Eigen::VectorX<T>& v, size_t block_idx,
                const std::vector<T>& values)
{
    size_t start = block_idx * 6;
    for (size_t i = 0; i < 6 && i < values.size(); ++i)
    {
        v(start + i) = values[i];
    }
}

template <typename T>
void set_block7(Eigen::VectorX<T>& v, size_t block_idx,
                const std::vector<T>& values)
{
    size_t start = block_idx * 7;
    for (size_t i = 0; i < 7 && i < values.size(); ++i)
    {
        v(start + i) = values[i];
    }
}

template <typename T>
void set_block4(Eigen::VectorX<T>& v, size_t block_idx,
                const std::vector<T>& values)
{
    size_t start = block_idx * 4;
    for (size_t i = 0; i < 4 && i < values.size(); ++i)
    {
        v(start + i) = values[i];
    }
}

template <typename T>
void set_body_velocity(Eigen::VectorX<T>& u, size_t body_idx,
                       const std::vector<T>& values)
{
    size_t start = body_idx * 6;
    for (size_t i = 0; i < 6 && i < values.size(); ++i)
    {
        u(start + i) = values[i];
    }
}

//template <typename prox_solver_functor>
template <typename T>
void run_2_body_prox_solver_test_eigen(
    T dummy_val = T(0)
    /*prox_solver_functor const& prox_solver*/)
{

    // Create vectors with appropriate sizes
    Eigen::VectorX<T> q(14);     // 2 bodies * 7
    Eigen::VectorX<T> Wdth(12);  // 2 bodies * 6
    Eigen::VectorX<T> u(12);     // 2 bodies * 6
    Eigen::VectorX<T> lambda(4); // 1 contact * 4
    Eigen::VectorX<T> g(4);      // 1 contact * 4
    Eigen::VectorX<T> e(4);      // 1 contact * 4
    Eigen::VectorX<T> mu(4);     // 1 contact * 4

    // Create sparse matrices
    Eigen::SparseMatrix<T> J(
        4, 12);  // 4 rows (1 contact * 4), 12 cols (2 bodies * 6)
    Eigen::SparseMatrix<T> WJT(
        12, 4); // 12 rows (2 bodies * 6), 4 cols (1 contact * 4)

    // Filling the vectors
    set_block7(q, 0, {0, 0, 0, 0, 0, 0, -1.0});
    set_block7(q, 1, {0, 2.0, 0, 0, 0, 0, -1.0});

    set_block6(u, 0, {0, 0, 0, 0, 0, 0});
    set_block6(u, 1, {0, 0, 0, 0, 0, 0});

    set_block6(Wdth, 0, {0, 0, 0, 0, 0, 0});
    set_block6(Wdth, 1, {0, -0.0981, 0, 0, 0, 0});

    set_block4(g, 0, {0, 0, 0, 0});
    set_block4(e, 0, {1, 0, 0, 0});
    set_block4(mu, 0, {0, 0.5, 0.5, 0.5});

    // Filling the matrices using triplets
    std::vector<Eigen::Triplet<T>> J_triplets;
    std::vector<Eigen::Triplet<T>> WJT_triplets;

    // J(0, 0) - body 0
    J_triplets.emplace_back(0, 0, 0);
    J_triplets.emplace_back(0, 1, -1);
    J_triplets.emplace_back(0, 2, 0);
    J_triplets.emplace_back(0, 3, 0);
    J_triplets.emplace_back(0, 4, 0);
    J_triplets.emplace_back(0, 5, 0);
    J_triplets.emplace_back(1, 0, -1);
    J_triplets.emplace_back(1, 1, 0);
    J_triplets.emplace_back(1, 2, 0);
    J_triplets.emplace_back(1, 3, 0);
    J_triplets.emplace_back(1, 4, 0);
    J_triplets.emplace_back(1, 5, 1);
    J_triplets.emplace_back(2, 0, 0);
    J_triplets.emplace_back(2, 1, 0);
    J_triplets.emplace_back(2, 2, -1);
    J_triplets.emplace_back(2, 3, -1);
    J_triplets.emplace_back(2, 4, 0);
    J_triplets.emplace_back(2, 5, 0);
    J_triplets.emplace_back(3, 0, 0);
    J_triplets.emplace_back(3, 1, 0);
    J_triplets.emplace_back(3, 2, 0);
    J_triplets.emplace_back(3, 3, 0);
    J_triplets.emplace_back(3, 4, -1);
    J_triplets.emplace_back(3, 5, 0);

    // J(0, 1) - body 1
    J_triplets.emplace_back(0, 6, 0);
    J_triplets.emplace_back(0, 7, 1);
    J_triplets.emplace_back(0, 8, 0);
    J_triplets.emplace_back(0, 9, 0);
    J_triplets.emplace_back(0, 10, 0);
    J_triplets.emplace_back(0, 11, 0);
    J_triplets.emplace_back(1, 6, 1);
    J_triplets.emplace_back(1, 7, 0);
    J_triplets.emplace_back(1, 8, 0);
    J_triplets.emplace_back(1, 9, 0);
    J_triplets.emplace_back(1, 10, 0);
    J_triplets.emplace_back(1, 11, 1);
    J_triplets.emplace_back(2, 6, 0);
    J_triplets.emplace_back(2, 7, 0);
    J_triplets.emplace_back(2, 8, 1);
    J_triplets.emplace_back(2, 9, -1);
    J_triplets.emplace_back(2, 10, 0);
    J_triplets.emplace_back(2, 11, 0);
    J_triplets.emplace_back(3, 6, 0);
    J_triplets.emplace_back(3, 7, 0);
    J_triplets.emplace_back(3, 8, 0);
    J_triplets.emplace_back(3, 9, 0);
    J_triplets.emplace_back(3, 10, 1);
    J_triplets.emplace_back(3, 11, 0);

    // WJT(0, 0) - body 0
    for (int i = 0; i < 24; ++i) { WJT_triplets.emplace_back(i % 6, i / 6, 0); }

    // WJT(1, 0) - body 1
    WJT_triplets.emplace_back(0, 0, 0);
    WJT_triplets.emplace_back(1, 0, 0.05);
    WJT_triplets.emplace_back(2, 0, 0);
    WJT_triplets.emplace_back(3, 0, 0);
    WJT_triplets.emplace_back(4, 0, 0.05);
    WJT_triplets.emplace_back(5, 0, 0);
    WJT_triplets.emplace_back(0, 1, 0);
    WJT_triplets.emplace_back(1, 1, 0);
    WJT_triplets.emplace_back(2, 1, 0);
    WJT_triplets.emplace_back(3, 1, 0);
    WJT_triplets.emplace_back(4, 1, 0);
    WJT_triplets.emplace_back(5, 1, 0.05);
    WJT_triplets.emplace_back(0, 2, 0);
    WJT_triplets.emplace_back(1, 2, 0);
    WJT_triplets.emplace_back(2, 2, 0);
    WJT_triplets.emplace_back(3, 2, -0.05);
    WJT_triplets.emplace_back(4, 2, 0);
    WJT_triplets.emplace_back(5, 2, 0);
    WJT_triplets.emplace_back(0, 3, 0);
    WJT_triplets.emplace_back(1, 3, 0.05);
    WJT_triplets.emplace_back(2, 3, 0);
    WJT_triplets.emplace_back(3, 3, 0.05);
    WJT_triplets.emplace_back(4, 3, 0);
    WJT_triplets.emplace_back(5, 3, 0);

    // Set the matrices from triplets
    J.setFromTriplets(J_triplets.begin(), J_triplets.end());
    WJT.setFromTriplets(WJT_triplets.begin(), WJT_triplets.end());
    T tolerance = 1e-1;
    BOOST_CHECK_CLOSE(lambda(0), 1.962f,
                      tolerance);  // Equivalent to lambda(0)[0]
    BOOST_CHECK_EQUAL(lambda(1), 0); // Equivalent to lambda(0)[1]
    BOOST_CHECK_EQUAL(lambda(2), 0); // Equivalent to lambda(0)[2]
    BOOST_CHECK_EQUAL(lambda(3), 0); // Equivalent to lambda(0)[3]

    // We now apply a small velocity in the direction of the x axis on the sphere.
    // This should result in friction. If the speed is 0.25, the impulse of the sphere
    // is P = m*v = 20*0.25 = 5. The bound on the friction impulse is
    // mu*l_n = 0.5*1.962 = 0.981.
    lambda.setZero(); // Equivalent to lambda.clear()
    set_body_velocity(u, 0, {0.25, 0, 0, 0, 0, 0}); // Set u(0)[0] = 0.25

    BOOST_CHECK(false);
    // prox_solver( ... fix me .... );

    BOOST_CHECK_CLOSE(lambda(0), 1.962f, tolerance);
    BOOST_CHECK_CLOSE(lambda(1), 0.981f, tolerance);
    BOOST_CHECK_EQUAL(lambda(2), 0);
    BOOST_CHECK_EQUAL(lambda(3), 0);

    // Applying a small velocity in the direction of the z axis on the sphere, should
    // result in friction. If the speed is 0.01, the impulse is P = m*v = 20*0.01 = 0.2.
    // This changes the bounds on the friction to (0.980, 0.0390)
    lambda.setZero();
    set_body_velocity(u, 0, {0.25, 0, 0.01, 0, 0, 0}); // Set u(0)[2] = 0.01

    BOOST_CHECK(false);
    // prox_solver( ... fix me .... );

    BOOST_CHECK_CLOSE(lambda(0), 1.962f, tolerance);
    BOOST_CHECK_CLOSE(lambda(1), 0.980f, tolerance);
    BOOST_CHECK_CLOSE(lambda(2), 0.0392f, tolerance);
    BOOST_CHECK_EQUAL(lambda(3), 0);

    // We now apply an angular velocity around the y axis on the sphere, to test Contensou friction
    // Setting w_y = 0.25, yields an upper bound of mu l_n = 0.981 (we assume the 2/3 factor to be
    // included in mu)
    lambda.setZero();
    set_body_velocity(u, 0, {0, 0, 0, 0, 0.25, 0}); // Set u(0)[4] = 0.25

    BOOST_CHECK(false);
    // prox_solver( ... fix me .... );

    BOOST_CHECK_CLOSE(lambda(0), 1.962f, tolerance);
    BOOST_CHECK_EQUAL(lambda(1), 0);
    BOOST_CHECK_EQUAL(lambda(2), 0);
    BOOST_CHECK_CLOSE(lambda(3), 0.981f, tolerance);
}

/**
 * This function build a simple test setup for use with a solver for the contact problem,
 * and performs a number of tests. It takes the math policy and the solver function as
 * template arguments.
 *
 * The setup has a unit sphere b_1 placed at (0,2,0), resting
 * on top of a 2x2x2 fixed box b_0 with its center at the origin. Another unit
 * sphere b_2 is present at (1.4142,3.4142,0) and is in contact with s_1. All
 * bodies are at rest.
 *
 * The sphere b_2 is expected to experience a normal velocity of magnitude
 * v = 0.5gt = 0.5*9.81*0.01 = 0.04905, with the impulse P = mv
 * = 20*0.04905 = 0.981. In the friction direction, the velocity is
 * expected to be similar, but the bound on the friction impulse is
 * mu*lambda_n = 0.5*0.981 = 0.4905.
 *
 * The sphere b_1 is affected by two contact points, with b_0 and b_1.
 * The speed of the sphere is g*t = -9.81*0.01 = -0.0981, and the impulse is
 * P = m*v = 20*-0.0981 = -1.962, so the magnitude of the normal impulse must
 * be close to 1.962.
 *
 * However, as the two contact points affect the same body, b_1, they will
 * affect each other. The upper sphere b_2 will push b_1 to the side, causing
 * a friction impulse between b_0 and b_1, smaller in magnitude than the friction
 * between b_1 and b_2, as it is this impulse that causes the movement of b_1.
 * As the lower sphere will move out of the way, the actual normal impulse
 * between b_1 and b_2 is likely to be smaller than computed above, however, the
 * normal impulse between b_0 and b_1 is likely to be larger than computed above.
 *
 * The solution is exected to be (3.428,0.5985,0,0) and (1.424,-0.6498,0,0)
 */
/*template <typename math_policy, typename prox_solver_functor>
void run_3_body_prox_solver_test(prox_solver_functor const& prox_solver)
{
    typename math_policy::vector7_type q(3);
    typename math_policy::vector6_type Wdth(3), u(3);
    typename math_policy::vector4_type lambda(2), g(2), e(2), mu(2);
    typename math_policy::compressed4x6_type J(2, 3, 4);
    typename math_policy::compressed6x4_type WJT(3, 2, 4);
    prox::SolverParams<math_policy> params;

  //filling the vectors
    q(0)[0] = 0;
    q(0)[1] = 0;
    q(0)[2] = 0;
    q(0)[3] = 0;
    q(0)[4] = 0;
    q(0)[5] = 0;
    q(0)[6] = -1.0;
    q(1)[0] = 0;
    q(1)[1] = 2.0;
    q(1)[2] = 0;
    q(1)[3] = 0;
    q(1)[4] = 0;
    q(1)[5] = 0;
    q(1)[6] = -1.0;
    q(2)[0] = 1.415;
    q(2)[1] = 3.415;
    q(2)[2] = 0;
    q(2)[3] = 0;
    q(1)[4] = 0;
    q(2)[5] = 0;
    q(2)[6] = -1.0;

    u(0)[0] = 0;
    u(0)[1] = 0;
    u(0)[2] = 0;
    u(0)[3] = 0;
    u(0)[4] = 0;
    u(0)[5] = 0;
    u(1)[0] = 0;
    u(1)[1] = 0;
    u(1)[2] = 0;
    u(1)[3] = 0;
    u(1)[4] = 0;
    u(1)[5] = 0;
    u(2)[0] = 0;
    u(2)[1] = 0;
    u(2)[2] = 0;
    u(2)[3] = 0;
    u(2)[4] = 0;
    u(2)[5] = 0;

    Wdth(0)[0] = 0;
    Wdth(0)[1] = 0;
    Wdth(0)[2] = 0;
    Wdth(0)[3] = 0;
    Wdth(0)[4] = 0;
    Wdth(0)[5] = 0;
    Wdth(1)[0] = 0;
    Wdth(1)[1] = -0.0981;
    Wdth(1)[2] = 0;
    Wdth(1)[3] = 0;
    Wdth(1)[4] = 0;
    Wdth(1)[5] = 0;
    Wdth(2)[0] = 0;
    Wdth(2)[1] = -0.0981;
    Wdth(1)[2] = 0;
    Wdth(2)[3] = 0;
    Wdth(2)[4] = 0;
    Wdth(2)[5] = 0;

    g(0)[0] = 0;
    g(0)[1] = 0;
    g(0)[2] = 0;
    g(0)[3] = 0;
    g(1)[0] = 0;
    g(1)[1] = 0;
    g(1)[2] = 0;
    g(1)[3] = 0;
    e(0)[0] = 1;
    e(0)[1] = 0;
    e(0)[2] = 0;
    e(0)[3] = 0;
    e(1)[0] = 1;
    e(1)[1] = 0;
    e(1)[2] = 0;
    e(1)[3] = 0;
    mu(0)[0] = 0;
    mu(0)[1] = 0.5;
    mu(0)[2] = 0.5;
    mu(0)[3] = 0.5;
    mu(1)[0] = 0;
    mu(1)[1] = 0.5;
    mu(1)[2] = 0.5;
    mu(1)[3] = 0.5;

  //filling the matrices
  //math_policy::block4x6_type & Jblock = J(0,0);
    J(0, 0)[0] = 0;
    J(0, 0)[1] = -1;
    J(0, 0)[2] = 0;
    J(0, 0)[3] = 0;
    J(0, 0)[4] = 0;
    J(0, 0)[5] = 0;
    J(0, 0)[6] = -1;
    J(0, 0)[7] = 0;
    J(0, 0)[8] = 0;
    J(0, 0)[9] = 0;
    J(0, 0)[10] = 0;
    J(0, 0)[11] = 1;
    J(0, 0)[12] = 0;
    J(0, 0)[13] = 0;
    J(0, 0)[14] = -1;
    J(0, 0)[15] = -1;
    J(0, 0)[16] = 0;
    J(0, 0)[17] = 0;
    J(0, 0)[18] = 0;
    J(0, 0)[19] = 0;
    J(0, 0)[20] = 0;
    J(0, 0)[21] = 0;
    J(0, 0)[22] = -1;
    J(0, 0)[23] = 0;

  //math_policy::block4x6_type & Jblock2 = J(0,1);
    J(0, 1)[0] = 0;
    J(0, 1)[1] = 1;
    J(0, 1)[2] = 0;
    J(0, 1)[3] = 0;
    J(0, 1)[4] = 0;
    J(0, 1)[5] = 0;
    J(0, 1)[6] = 1;
    J(0, 1)[7] = 0;
    J(0, 1)[8] = 0;
    J(0, 1)[9] = 0;
    J(0, 1)[10] = 0;
    J(0, 1)[11] = 1;
    J(0, 1)[12] = 0;
    J(0, 1)[13] = 0;
    J(0, 1)[14] = 1;
    J(0, 1)[15] = -1;
    J(0, 1)[16] = 0;
    J(0, 1)[17] = 0;
    J(0, 1)[18] = 0;
    J(0, 1)[19] = 0;
    J(0, 1)[20] = 0;
    J(0, 1)[21] = 0;
    J(0, 1)[22] = 1;
    J(0, 1)[23] = 0;

  //math_policy::block4x6_type & Jblock = J(1,1);
  //math_policy::real_type c0 = sqrt(2.0f);
    typename math_policy::real_type c1 = 1.0f / sqrt(2.0f);
    J(1, 1)[0] = -c1;
    J(1, 1)[1] = -c1;
    J(1, 1)[2] = 0;
    J(1, 1)[3] = 0;
    J(1, 1)[4] = 0;
    J(1, 1)[5] = 0;
    J(1, 1)[6] = -c1;
    J(1, 1)[7] = c1;
    J(1, 1)[8] = 0;
    J(1, 1)[9] = 0;
    J(1, 1)[10] = 0;
    J(1, 1)[11] = 1;
    J(1, 1)[12] = 0;
    J(1, 1)[13] = 0;
    J(1, 1)[14] = -1;
    J(1, 1)[15] = -c1;
    J(1, 1)[16] = c1;
    J(1, 1)[17] = 0;
    J(1, 1)[18] = 0;
    J(1, 1)[19] = 0;
    J(1, 1)[20] = 0;
    J(1, 1)[21] = -c1;
    J(1, 1)[22] = -c1;
    J(1, 1)[23] = 0;

  //math_policy::block4x6_type & Jblock2 = J(1,2);
    J(1, 2)[0] = c1;
    J(1, 2)[1] = c1;
    J(1, 2)[2] = 0;
    J(1, 2)[3] = 0;
    J(1, 2)[4] = 0;
    J(1, 2)[5] = 0;
    J(1, 2)[6] = c1;
    J(1, 2)[7] = -c1;
    J(1, 2)[8] = 0;
    J(1, 2)[9] = 0;
    J(1, 2)[10] = 0;
    J(1, 2)[11] = 1;
    J(1, 2)[12] = 0;
    J(1, 2)[13] = 0;
    J(1, 2)[14] = 1;
    J(1, 2)[15] = -c1;
    J(1, 2)[16] = c1;
    J(1, 2)[17] = 0;
    J(1, 2)[18] = 0;
    J(1, 2)[19] = 0;
    J(1, 2)[20] = 0;
    J(1, 2)[21] = c1;
    J(1, 2)[22] = c1;
    J(1, 2)[23] = 0;

  //math_policy::block6x4_type & WJTblock = WJT(0,0);
    WJT(0, 0)[0] = 0;
    WJT(0, 0)[1] = 0;
    WJT(0, 0)[2] = 0;
    WJT(0, 0)[3] = 0;
    WJT(0, 0)[4] = 0;
    WJT(0, 0)[5] = 0;
    WJT(0, 0)[6] = 0;
    WJT(0, 0)[7] = 0;
    WJT(0, 0)[8] = 0;
    WJT(0, 0)[9] = 0;
    WJT(0, 0)[10] = 0;
    WJT(0, 0)[11] = 0;
    WJT(0, 0)[12] = 0;
    WJT(0, 0)[13] = 0;
    WJT(0, 0)[14] = 0;
    WJT(0, 0)[15] = 0;
    WJT(0, 0)[16] = 0;
    WJT(0, 0)[17] = 0;
    WJT(0, 0)[18] = 0;
    WJT(0, 0)[19] = 0;
    WJT(0, 0)[20] = 0;
    WJT(0, 0)[21] = 0;
    WJT(0, 0)[22] = 0;
    WJT(0, 0)[23] = 0;
  //WJTblock = WJT(1,0);
    WJT(1, 0)[0] = 0;
    WJT(1, 0)[1] = 0.05;
    WJT(1, 0)[2] = 0;
    WJT(1, 0)[3] = 0;
    WJT(1, 0)[4] = 0.05;
    WJT(1, 0)[5] = 0;
    WJT(1, 0)[6] = 0;
    WJT(1, 0)[7] = 0;
    WJT(1, 0)[8] = 0;
    WJT(1, 0)[9] = 0;
    WJT(1, 0)[10] = 0.05;
    WJT(1, 0)[11] = 0;
    WJT(1, 0)[12] = 0;
    WJT(1, 0)[13] = 0;
    WJT(1, 0)[14] = -0.05;
    WJT(1, 0)[15] = 0;
    WJT(1, 0)[16] = 0;
    WJT(1, 0)[17] = 0;
    WJT(1, 0)[18] = 0;
    WJT(1, 0)[19] = 0.05;
    WJT(1, 0)[20] = 0;
    WJT(1, 0)[21] = 0.05;
    WJT(1, 0)[22] = 0;
    WJT(1, 0)[23] = 0;
  //math_policy::block6x4_type & WJTblock = WJT(1,1);
  //math_policy::real_type c2 = c0/20.0f;
    typename math_policy::real_type c3 = c1 / 20.0f;
    WJT(1, 1)[0] = -c3;
    WJT(1, 1)[1] = -c3;
    WJT(1, 1)[2] = 0;
    WJT(1, 1)[3] = 0;
    WJT(1, 1)[4] = -c3;
    WJT(1, 1)[5] = c3;
    WJT(1, 1)[6] = 0;
    WJT(1, 1)[7] = 0;
    WJT(1, 1)[8] = 0;
    WJT(1, 1)[9] = 0;
    WJT(1, 1)[10] = -0.05;
    WJT(1, 1)[11] = 0;
    WJT(1, 1)[12] = 0;
    WJT(1, 1)[13] = 0;
    WJT(1, 1)[14] = -c3;
    WJT(1, 1)[15] = -c3;
    WJT(1, 1)[16] = 0;
    WJT(1, 1)[17] = 0;
    WJT(1, 1)[18] = c3;
    WJT(1, 1)[19] = -c3;
    WJT(1, 1)[20] = 0;
    WJT(1, 1)[21] = 0.05;
    WJT(1, 1)[22] = 0;
    WJT(1, 1)[23] = 0;
  //WJTblock = WJT(2,1);
    WJT(2, 1)[0] = c3;
    WJT(2, 1)[1] = c3;
    WJT(2, 1)[2] = 0;
    WJT(2, 1)[3] = 0;
    WJT(2, 1)[4] = c3;
    WJT(2, 1)[5] = -c3;
    WJT(2, 1)[6] = 0;
    WJT(2, 1)[7] = 0;
    WJT(2, 1)[8] = 0;
    WJT(2, 1)[9] = 0;
    WJT(2, 1)[10] = 0.05;
    WJT(2, 1)[11] = 0;
    WJT(2, 1)[12] = 0;
    WJT(2, 1)[13] = 0;
    WJT(2, 1)[14] = c3;
    WJT(2, 1)[15] = c3;
    WJT(2, 1)[16] = 0;
    WJT(2, 1)[17] = 0;
    WJT(2, 1)[18] = c3;
    WJT(2, 1)[19] = c3;
    WJT(2, 1)[20] = 0;
    WJT(2, 1)[21] = 0.05;
    WJT(2, 1)[22] = 0;
    WJT(2, 1)[23] = 0;

    BOOST_CHECK(false);
  //  prox_solver( ... fix me .... );

    typename math_policy::real_type tolerance = 1e-1; // that is 0.1 percent, an error on the fourth decimal
    BOOST_CHECK_CLOSE(lambda(0)[0], 3.428f, tolerance);
    BOOST_CHECK_CLOSE(lambda(0)[1], 0.5985f, tolerance);
    BOOST_CHECK_EQUAL(lambda(0)[2], 0);
    BOOST_CHECK_EQUAL(lambda(0)[3], 0);

    BOOST_CHECK_CLOSE(lambda(1)[0], 1.424f, tolerance);
    BOOST_CHECK_CLOSE(lambda(1)[1], -0.6498f, tolerance);
    BOOST_CHECK_EQUAL(lambda(1)[2], 0);
    BOOST_CHECK_EQUAL(lambda(1)[3], 0);
}*/

template <typename T>
void run_3_body_prox_solver_test_eigen(T i_am_a_dummy = T(0))
{

    // Create vectors with appropriate sizes
    Eigen::VectorX<T> q(21); // 3 bodies * 7
    Eigen::VectorX<T> Wdth(18); // 3 bodies * 6
    Eigen::VectorX<T> u(18); // 3 bodies * 6
    Eigen::VectorX<T> lambda(8); // 2 contacts * 4
    Eigen::VectorX<T> g(8); // 2 contacts * 4
    Eigen::VectorX<T> e(8); // 2 contacts * 4
    Eigen::VectorX<T> mu(8); // 2 contacts * 4

    // Create sparse matrices
    Eigen::SparseMatrix<T> J(
        8, 18); // 8 rows (2 contacts * 4), 18 cols (3 bodies * 6)
    Eigen::SparseMatrix<T> WJT(
        18, 8); // 18 rows (3 bodies * 6), 8 cols (2 contacts * 4)

    // Solver parameters
    prox::SolverParams<T> params;

    // Filling the vectors
    set_block7(q, 0, {0, 0, 0, 0, 0, 0, -1.0});
    set_block7(q, 1, {0, 2.0, 0, 0, 0, 0, -1.0});
    set_block7(q, 2, {1.415, 3.415, 0, 0, 0, 0, -1.0});

    set_block6(u, 0, {0, 0, 0, 0, 0, 0});
    set_block6(u, 1, {0, 0, 0, 0, 0, 0});
    set_block6(u, 2, {0, 0, 0, 0, 0, 0});

    set_block6(Wdth, 0, {0, 0, 0, 0, 0, 0});
    set_block6(Wdth, 1, {0, -0.0981, 0, 0, 0, 0});
    set_block6(Wdth, 2, {0, -0.0981, 0, 0, 0, 0});

    set_block4(g, 0, {0, 0, 0, 0});
    set_block4(g, 1, {0, 0, 0, 0});
    set_block4(e, 0, {1, 0, 0, 0});
    set_block4(e, 1, {1, 0, 0, 0});
    set_block4(mu, 0, {0, 0.5, 0.5, 0.5});
    set_block4(mu, 1, {0, 0.5, 0.5, 0.5});

    // Filling the matrices using triplets
    std::vector<Eigen::Triplet<T>> J_triplets;
    std::vector<Eigen::Triplet<T>> WJT_triplets;

    // J(0, 0) - contact 0, body 0
    J_triplets.emplace_back(0, 0, 0);
    J_triplets.emplace_back(0, 1, -1);
    J_triplets.emplace_back(0, 2, 0);
    J_triplets.emplace_back(0, 3, 0);
    J_triplets.emplace_back(0, 4, 0);
    J_triplets.emplace_back(0, 5, 0);
    J_triplets.emplace_back(1, 0, -1);
    J_triplets.emplace_back(1, 1, 0);
    J_triplets.emplace_back(1, 2, 0);
    J_triplets.emplace_back(1, 3, 0);
    J_triplets.emplace_back(1, 4, 0);
    J_triplets.emplace_back(1, 5, 1);
    J_triplets.emplace_back(2, 0, 0);
    J_triplets.emplace_back(2, 1, 0);
    J_triplets.emplace_back(2, 2, -1);
    J_triplets.emplace_back(2, 3, -1);
    J_triplets.emplace_back(2, 4, 0);
    J_triplets.emplace_back(2, 5, 0);
    J_triplets.emplace_back(3, 0, 0);
    J_triplets.emplace_back(3, 1, 0);
    J_triplets.emplace_back(3, 2, 0);
    J_triplets.emplace_back(3, 3, 0);
    J_triplets.emplace_back(3, 4, -1);
    J_triplets.emplace_back(3, 5, 0);

    // J(0, 1) - contact 0, body 1
    J_triplets.emplace_back(0, 6, 0);
    J_triplets.emplace_back(0, 7, 1);
    J_triplets.emplace_back(0, 8, 0);
    J_triplets.emplace_back(0, 9, 0);
    J_triplets.emplace_back(0, 10, 0);
    J_triplets.emplace_back(0, 11, 0);
    J_triplets.emplace_back(1, 6, 1);
    J_triplets.emplace_back(1, 7, 0);
    J_triplets.emplace_back(1, 8, 0);
    J_triplets.emplace_back(1, 9, 0);
    J_triplets.emplace_back(1, 10, 0);
    J_triplets.emplace_back(1, 11, 1);
    J_triplets.emplace_back(2, 6, 0);
    J_triplets.emplace_back(2, 7, 0);
    J_triplets.emplace_back(2, 8, 1);
    J_triplets.emplace_back(2, 9, -1);
    J_triplets.emplace_back(2, 10, 0);
    J_triplets.emplace_back(2, 11, 0);
    J_triplets.emplace_back(3, 6, 0);
    J_triplets.emplace_back(3, 7, 0);
    J_triplets.emplace_back(3, 8, 0);
    J_triplets.emplace_back(3, 9, 0);
    J_triplets.emplace_back(3, 10, 1);
    J_triplets.emplace_back(3, 11, 0);

    // J(1, 1) - contact 1, body 1
    T c1 = 1.0f / std::sqrt(2.0f);
    J_triplets.emplace_back(4, 6, -c1);
    J_triplets.emplace_back(4, 7, -c1);
    J_triplets.emplace_back(4, 8, 0);
    J_triplets.emplace_back(4, 9, 0);
    J_triplets.emplace_back(4, 10, 0);
    J_triplets.emplace_back(4, 11, 0);
    J_triplets.emplace_back(5, 6, -c1);
    J_triplets.emplace_back(5, 7, c1);
    J_triplets.emplace_back(5, 8, 0);
    J_triplets.emplace_back(5, 9, 0);
    J_triplets.emplace_back(5, 10, 0);
    J_triplets.emplace_back(5, 11, 1);
    J_triplets.emplace_back(6, 6, 0);
    J_triplets.emplace_back(6, 7, 0);
    J_triplets.emplace_back(6, 8, -1);
    J_triplets.emplace_back(6, 9, -c1);
    J_triplets.emplace_back(6, 10, c1);
    J_triplets.emplace_back(6, 11, 0);
    J_triplets.emplace_back(7, 6, 0);
    J_triplets.emplace_back(7, 7, 0);
    J_triplets.emplace_back(7, 8, 0);
    J_triplets.emplace_back(7, 9, 0);
    J_triplets.emplace_back(7, 10, 0);
    J_triplets.emplace_back(7, 11, -c1);
    J_triplets.emplace_back(
        7, 6, -c1); // Note: This might be a duplicate, verify original indices

    // J(1, 2) - contact 1, body 2
    J_triplets.emplace_back(4, 12, c1);
    J_triplets.emplace_back(4, 13, c1);
    J_triplets.emplace_back(4, 14, 0);
    J_triplets.emplace_back(4, 15, 0);
    J_triplets.emplace_back(4, 16, 0);
    J_triplets.emplace_back(4, 17, 0);
    J_triplets.emplace_back(5, 12, c1);
    J_triplets.emplace_back(5, 13, -c1);
    J_triplets.emplace_back(5, 14, 0);
    J_triplets.emplace_back(5, 15, 0);
    J_triplets.emplace_back(5, 16, 0);
    J_triplets.emplace_back(5, 17, 1);
    J_triplets.emplace_back(6, 12, 0);
    J_triplets.emplace_back(6, 13, 0);
    J_triplets.emplace_back(6, 14, 1);
    J_triplets.emplace_back(6, 15, -c1);
    J_triplets.emplace_back(6, 16, c1);
    J_triplets.emplace_back(6, 17, 0);
    J_triplets.emplace_back(7, 12, 0);
    J_triplets.emplace_back(7, 13, 0);
    J_triplets.emplace_back(7, 14, 0);
    J_triplets.emplace_back(7, 15, 0);
    J_triplets.emplace_back(7, 16, 0);
    J_triplets.emplace_back(7, 17, c1);
    J_triplets.emplace_back(
        7, 12, c1); // Note: This might be a duplicate, verify original indices

    // WJT(0, 0) - body 0, contact 0
    for (int i = 0; i < 24; ++i) { WJT_triplets.emplace_back(i % 6, i / 6, 0); }

    // WJT(1, 0) - body 1, contact 0
    WJT_triplets.emplace_back(0, 0, 0);
    WJT_triplets.emplace_back(1, 0, 0.05);
    WJT_triplets.emplace_back(2, 0, 0);
    WJT_triplets.emplace_back(3, 0, 0);
    WJT_triplets.emplace_back(4, 0, 0.05);
    WJT_triplets.emplace_back(5, 0, 0);
    WJT_triplets.emplace_back(0, 1, 0);
    WJT_triplets.emplace_back(1, 1, 0);
    WJT_triplets.emplace_back(2, 1, 0);
    WJT_triplets.emplace_back(3, 1, 0);
    WJT_triplets.emplace_back(4, 1, 0);
    WJT_triplets.emplace_back(5, 1, 0.05);
    WJT_triplets.emplace_back(0, 2, 0);
    WJT_triplets.emplace_back(1, 2, 0);
    WJT_triplets.emplace_back(2, 2, 0);
    WJT_triplets.emplace_back(3, 2, -0.05);
    WJT_triplets.emplace_back(4, 2, 0);
    WJT_triplets.emplace_back(5, 2, 0);
    WJT_triplets.emplace_back(0, 3, 0);
    WJT_triplets.emplace_back(1, 3, 0.05);
    WJT_triplets.emplace_back(2, 3, 0);
    WJT_triplets.emplace_back(3, 3, 0.05);
    WJT_triplets.emplace_back(4, 3, 0);
    WJT_triplets.emplace_back(5, 3, 0);

    // WJT(1, 1) - body 1, contact 1
    T c3 = c1 / 20.0f;
    WJT_triplets.emplace_back(6, 4, -c3);
    WJT_triplets.emplace_back(7, 4, -c3);
    WJT_triplets.emplace_back(8, 4, 0);
    WJT_triplets.emplace_back(9, 4, 0);
    WJT_triplets.emplace_back(10, 4, -c3);
    WJT_triplets.emplace_back(11, 4, c3);
    WJT_triplets.emplace_back(6, 5, 0);
    WJT_triplets.emplace_back(7, 5, 0);
    WJT_triplets.emplace_back(8, 5, 0);
    WJT_triplets.emplace_back(9, 5, 0);
    WJT_triplets.emplace_back(10, 5, -0.05);
    WJT_triplets.emplace_back(11, 5, 0);
    WJT_triplets.emplace_back(6, 6, 0);
    WJT_triplets.emplace_back(7, 6, 0);
    WJT_triplets.emplace_back(8, 6, 0);
    WJT_triplets.emplace_back(9, 6, 0);
    WJT_triplets.emplace_back(10, 6, -c3);
    WJT_triplets.emplace_back(11, 6, -c3);
    WJT_triplets.emplace_back(6, 7, 0);
    WJT_triplets.emplace_back(7, 7, 0);
    WJT_triplets.emplace_back(8, 7, c3);
    WJT_triplets.emplace_back(9, 7, -c3);
    WJT_triplets.emplace_back(10, 7, 0);
    WJT_triplets.emplace_back(11, 7, 0.05);
    WJT_triplets.emplace_back(11, 7, 0); // Note: Verify original indices

    // WJT(2, 1) - body 2, contact 1
    WJT_triplets.emplace_back(12, 4, c3);
    WJT_triplets.emplace_back(13, 4, c3);
    WJT_triplets.emplace_back(14, 4, 0);
    WJT_triplets.emplace_back(15, 4, 0);
    WJT_triplets.emplace_back(16, 4, c3);
    WJT_triplets.emplace_back(17, 4, -c3);
    WJT_triplets.emplace_back(12, 5, 0);
    WJT_triplets.emplace_back(13, 5, 0);
    WJT_triplets.emplace_back(14, 5, 0);
    WJT_triplets.emplace_back(15, 5, 0);
    WJT_triplets.emplace_back(16, 5, 0.05);
    WJT_triplets.emplace_back(17, 5, 0);
    WJT_triplets.emplace_back(12, 6, 0);
    WJT_triplets.emplace_back(13, 6, 0);
    WJT_triplets.emplace_back(14, 6, 0);
    WJT_triplets.emplace_back(15, 6, 0);
    WJT_triplets.emplace_back(16, 6, c3);
    WJT_triplets.emplace_back(17, 6, c3);
    WJT_triplets.emplace_back(12, 7, 0);
    WJT_triplets.emplace_back(13, 7, 0);
    WJT_triplets.emplace_back(14, 7, c3);
    WJT_triplets.emplace_back(15, 7, c3);
    WJT_triplets.emplace_back(16, 7, 0);
    WJT_triplets.emplace_back(17, 7, 0.05);
    WJT_triplets.emplace_back(17, 7, 0); // Note: Verify original indices

    // Set the matrices from triplets
    J.setFromTriplets(J_triplets.begin(), J_triplets.end());
    WJT.setFromTriplets(WJT_triplets.begin(), WJT_triplets.end());

    BOOST_CHECK(false);
    // prox_solver(J, WJT, u, e, g, mu, lambda, params);

    T tolerance = 1e-1; // that is 0.1 percent, an error on the fourth decimal
    BOOST_CHECK_CLOSE(lambda(0), 3.428f, tolerance); // lambda(0)[0]
    BOOST_CHECK_CLOSE(lambda(1), 0.5985f, tolerance); // lambda(0)[1]
    BOOST_CHECK_EQUAL(lambda(2), 0); // lambda(0)[2]
    BOOST_CHECK_EQUAL(lambda(3), 0); // lambda(0)[3]

    BOOST_CHECK_CLOSE(lambda(4), 1.424f, tolerance); // lambda(1)[0]
    BOOST_CHECK_CLOSE(lambda(5), -0.6498f, tolerance); // lambda(1)[1]
    BOOST_CHECK_EQUAL(lambda(6), 0); // lambda(1)[2]
    BOOST_CHECK_EQUAL(lambda(7), 0); // lambda(1)[3]
}

/**
 * This function build a 4 body test setup for use with a solver for the contact problem,
 * and performs a number of tests. It takes the math policy as template arguments, and
 * the solver and number of cores to use as function arguments.
 *
 * The setup has a unit sphere b_1 placed at (0,2,0), resting
 * on top of a 2x2x2 fixed box b_0 with its center at the origin. Another unit
 * sphere b_2 is present at (1.4142,3.4142,0) and is in contact with b_1. A 3rd
 * unit sphere is placed at ( -0.517638,3.93185,0.0), which places it in contact
 * with b_1 and b_2. All bodies are at rest.
 *
 * The sphere b_2 is expected to experience a normal velocity of magnitude
 * v = 0.5gt = 0.5*9.81*0.01 = 0.04905, with the impulse P = mv
 * = 20*0.04905 = 0.981. In the friction direction, the velocity is
 * expected to be similar, but the bound on the friction impulse is
 * mu*lambda_n = 0.5*0.981 = 0.4905.
 *
 * The sphere b_1 is affected by two contact points, with b_0 and b_1.
 * The speed of the sphere is g*t = -9.81*0.01 = -0.0981, and the impulse is
 * P = m*v = 20*-0.0981 = -1.962, so the magnitude of the normal impulse must
 * be close to 1.962.
 *
 * However, as the two contact points affect the same body, b_1, they will
 * affect each other. The upper sphere b_2 will push b_1 to the side, causing
 * a friction impulse between b_0 and b_1, equal in magnitude to the friction
 * between b_1 and b_2, as it is this impulse that causes the movement of b_1.
 * As the lower sphere will move out of the way, the actual normal impulse
 * between b_1 and b_2 is likely to be smaller than computed above.
 *
 * The 3rd sphere is expected to push b_2 in the x direction and add normal
 * impulse to b_1.
 *
 * The solution is exected to be (5.327,0.3471,0,0), (1.410,-0.6658,0,0),
 *(1.887,0.2860,0,0) and (0,0,0,0)
 */
/*
template <typename math_policy, typename prox_solver_functor>
void run_4_body_prox_solver_test(prox_solver_functor const& prox_solver)
{
    typename math_policy::vector7_type q(4);
    typename math_policy::vector6_type Wdth(4), u(4);
    typename math_policy::vector4_type lambda(4), g(4), e(4), mu(4);
    typename math_policy::compressed4x6_type J(4, 4, 8);
    typename math_policy::compressed6x4_type WJT(4, 4, 8);
    prox::SolverParams<math_policy> params;

  //filling the vectors
    q(0)[0] = 0;
    q(0)[1] = 0;
    q(0)[2] = 0;
    q(0)[3] = 1.0;
    q(0)[4] = 0;
    q(0)[5] = 0;
    q(0)[6] = 0;
    q(1)[0] = 0;
    q(1)[1] = 2.0;
    q(1)[2] = 0;
    q(1)[3] = 1.0;
    q(1)[4] = 0;
    q(1)[5] = 0;
    q(1)[6] = 0;
    q(2)[0] = 1.4142;
    q(2)[1] = 3.4142;
    q(1)[2] = 0;
    q(2)[3] = 1.0;
    q(2)[4] = 0;
    q(2)[5] = 0;
    q(2)[6] = 0;
    q(3)[0] = -0.51764;
    q(3)[1] = 3.9319;
    q(3)[2] = 0;
    q(3)[3] = 1.0;
    q(3)[4] = 0;
    q(3)[5] = 0;
    q(3)[6] = 0;

    u(0)[0] = 0;
    u(0)[1] = 0;
    u(0)[2] = 0;
    u(0)[3] = 0;
    u(0)[4] = 0;
    u(0)[5] = 0;
    u(1)[0] = 0;
    u(1)[1] = 0;
    u(1)[2] = 0;
    u(1)[3] = 0;
    u(1)[4] = 0;
    u(1)[5] = 0;
    u(2)[0] = 0;
    u(2)[1] = 0;
    u(2)[2] = 0;
    u(2)[3] = 0;
    u(2)[4] = 0;
    u(2)[5] = 0;
    u(3)[0] = 0;
    u(3)[1] = 0;
    u(3)[2] = 0;
    u(3)[3] = 0;
    u(3)[4] = 0;
    u(3)[5] = 0;

    Wdth(0)[0] = 0;
    Wdth(0)[1] = 0;
    Wdth(0)[2] = 0;
    Wdth(0)[3] = 0;
    Wdth(0)[4] = 0;
    Wdth(0)[5] = 0;
    Wdth(1)[0] = 0;
    Wdth(1)[1] = -0.0981;
    Wdth(1)[2] = 0;
    Wdth(1)[3] = 0;
    Wdth(1)[4] = 0;
    Wdth(1)[5] = 0;
    Wdth(2)[0] = 0;
    Wdth(2)[1] = -0.0981;
    Wdth(1)[2] = 0;
    Wdth(2)[3] = 0;
    Wdth(2)[4] = 0;
    Wdth(2)[5] = 0;
    Wdth(3)[0] = 0;
    Wdth(3)[1] = -0.0981;
    Wdth(3)[2] = 0;
    Wdth(3)[3] = 0;
    Wdth(3)[4] = 0;
    Wdth(3)[5] = 0;

    g(0)[0] = 0;
    g(0)[1] = 0;
    g(0)[2] = 0;
    g(0)[3] = 0;
    g(1)[0] = 0;
    g(1)[1] = 0;
    g(1)[2] = 0;
    g(1)[3] = 0;
    g(2)[0] = 0;
    g(2)[1] = 0;
    g(2)[2] = 0;
    g(2)[3] = 0;
    g(3)[0] = 0;
    g(3)[1] = 0;
    g(3)[2] = 0;
    g(3)[3] = 0;

    e(0)[0] = 1;
    e(0)[1] = 0;
    e(0)[2] = 0;
    e(0)[3] = 0;
    e(1)[0] = 1;
    e(1)[1] = 0;
    e(1)[2] = 0;
    e(1)[3] = 0;
    e(2)[0] = 1;
    e(2)[1] = 0;
    e(2)[2] = 0;
    e(2)[3] = 0;
    e(3)[0] = 1;
    e(3)[1] = 0;
    e(3)[2] = 0;
    e(3)[3] = 0;

    mu(0)[0] = 0;
    mu(0)[1] = 0.5;
    mu(0)[2] = 0.5;
    mu(0)[3] = 0.5;
    mu(1)[0] = 0;
    mu(1)[1] = 0.5;
    mu(1)[2] = 0.5;
    mu(1)[3] = 0.5;
    mu(2)[0] = 0;
    mu(2)[1] = 0.5;
    mu(2)[2] = 0.5;
    mu(2)[3] = 0.5;
    mu(3)[0] = 0;
    mu(3)[1] = 0.5;
    mu(3)[2] = 0.5;
    mu(3)[3] = 0.5;

  //filling the matrices
  //math_policy::block4x6_type & Jblock = J(0,0);
    J(0, 0)[0] = 0;
    J(0, 0)[1] = -1;
    J(0, 0)[2] = 0;
    J(0, 0)[3] = 0;
    J(0, 0)[4] = 0;
    J(0, 0)[5] = 0;
    J(0, 0)[6] = -1;
    J(0, 0)[7] = 0;
    J(0, 0)[8] = 0;
    J(0, 0)[9] = 0;
    J(0, 0)[10] = 0;
    J(0, 0)[11] = 1;
    J(0, 0)[12] = 0;
    J(0, 0)[13] = 0;
    J(0, 0)[14] = -1;
    J(0, 0)[15] = -1;
    J(0, 0)[16] = 0;
    J(0, 0)[17] = 0;
    J(0, 0)[18] = 0;
    J(0, 0)[19] = 0;
    J(0, 0)[20] = 0;
    J(0, 0)[21] = 0;
    J(0, 0)[22] = -1;
    J(0, 0)[23] = 0;

  //math_policy::block4x6_type & Jblock2 = J(0,1);
    J(0, 1)[0] = 0;
    J(0, 1)[1] = 1;
    J(0, 1)[2] = 0;
    J(0, 1)[3] = 0;
    J(0, 1)[4] = 0;
    J(0, 1)[5] = 0;
    J(0, 1)[6] = 1;
    J(0, 1)[7] = 0;
    J(0, 1)[8] = 0;
    J(0, 1)[9] = 0;
    J(0, 1)[10] = 0;
    J(0, 1)[11] = 1;
    J(0, 1)[12] = 0;
    J(0, 1)[13] = 0;
    J(0, 1)[14] = 1;
    J(0, 1)[15] = -1;
    J(0, 1)[16] = 0;
    J(0, 1)[17] = 0;
    J(0, 1)[18] = 0;
    J(0, 1)[19] = 0;
    J(0, 1)[20] = 0;
    J(0, 1)[21] = 0;
    J(0, 1)[22] = 1;
    J(0, 1)[23] = 0;

  //math_policy::block4x6_type & Jblock = J(1,1);
  //math_policy::real_type c0 = sqrt(2.0f);
    typename math_policy::real_type c1 = 1.0f / sqrt(2.0f);
    J(1, 1)[0] = -c1;
    J(1, 1)[1] = -c1;
    J(1, 1)[2] = 0;
    J(1, 1)[3] = 0;
    J(1, 1)[4] = 0;
    J(1, 1)[5] = 0;
    J(1, 1)[6] = -c1;
    J(1, 1)[7] = c1;
    J(1, 1)[8] = 0;
    J(1, 1)[9] = 0;
    J(1, 1)[10] = 0;
    J(1, 1)[11] = 1;
    J(1, 1)[12] = 0;
    J(1, 1)[13] = 0;
    J(1, 1)[14] = -1;
    J(1, 1)[15] = -c1;
    J(1, 1)[16] = c1;
    J(1, 1)[17] = 0;
    J(1, 1)[18] = 0;
    J(1, 1)[19] = 0;
    J(1, 1)[20] = 0;
    J(1, 1)[21] = -c1;
    J(1, 1)[22] = -c1;
    J(1, 1)[23] = 0;

  //math_policy::block4x6_type & Jblock2 = J(1,2);
    J(1, 2)[0] = c1;
    J(1, 2)[1] = c1;
    J(1, 2)[2] = 0;
    J(1, 2)[3] = 0;
    J(1, 2)[4] = 0;
    J(1, 2)[5] = 0;
    J(1, 2)[6] = c1;
    J(1, 2)[7] = -c1;
    J(1, 2)[8] = 0;
    J(1, 2)[9] = 0;
    J(1, 2)[10] = 0;
    J(1, 2)[11] = 1;
    J(1, 2)[12] = 0;
    J(1, 2)[13] = 0;
    J(1, 2)[14] = 1;
    J(1, 2)[15] = -c1;
    J(1, 2)[16] = c1;
    J(1, 2)[17] = 0;
    J(1, 2)[18] = 0;
    J(1, 2)[19] = 0;
    J(1, 2)[20] = 0;
    J(1, 2)[21] = c1;
    J(1, 2)[22] = c1;
    J(1, 2)[23] = 0;

  // cp 3
    typename math_policy::real_type c30 = 0.258849;//0.259;
    typename math_policy::real_type c31 = 0.965927;//0.966;
    J(2, 1)[0] = c30;
    J(2, 1)[1] = -c31;
    J(2, 1)[2] = 0;
    J(2, 1)[3] = 0;
    J(2, 1)[4] = 0;
    J(2, 1)[5] = 0;
    J(2, 1)[6] = -c31;
    J(2, 1)[7] = -c30;
    J(2, 1)[8] = 0;
    J(2, 1)[9] = 0;
    J(2, 1)[10] = 0;
    J(2, 1)[11] = 1;
    J(2, 1)[12] = 0;
    J(2, 1)[13] = 0;
    J(2, 1)[14] = -1;
    J(2, 1)[15] = -c31;
    J(2, 1)[16] = -c30;
    J(2, 1)[17] = 0;
    J(2, 1)[18] = 0;
    J(2, 1)[19] = 0;
    J(2, 1)[20] = 0;
    J(2, 1)[21] = c30;
    J(2, 1)[22] = -c31;
    J(2, 1)[23] = 0;

  //math_policy::block4x6_type & Jblock2 = J(1,2);
    J(2, 3)[0] = -c30;
    J(2, 3)[1] = c31;
    J(2, 3)[2] = 0;
    J(2, 3)[3] = 0;
    J(2, 3)[4] = 0;
    J(2, 3)[5] = 0;
    J(2, 3)[6] = c31;
    J(2, 3)[7] = c30;
    J(2, 3)[8] = 0;
    J(2, 3)[9] = 0;
    J(2, 3)[10] = 0;
    J(2, 3)[11] = 1;
    J(2, 3)[12] = 0;
    J(2, 3)[13] = 0;
    J(2, 3)[14] = 1;
    J(2, 3)[15] = -c31;
    J(2, 3)[16] = -c30;
    J(2, 3)[17] = 0;
    J(2, 3)[18] = 0;
    J(2, 3)[19] = 0;
    J(2, 3)[20] = 0;
    J(2, 3)[21] = -c30;
    J(2, 3)[22] = c31;
    J(2, 3)[23] = 0;

  // cp4
    typename math_policy::real_type c40 = 0.258849;//0.259;
    typename math_policy::real_type c41 = 0.965927;//0.966;
    J(3, 2)[0] = c41;
    J(3, 2)[1] = -c40;
    J(3, 2)[2] = 0;
    J(3, 2)[3] = 0;
    J(3, 2)[4] = 0;
    J(3, 2)[5] = 0;
    J(3, 2)[6] = 0;
    J(3, 2)[7] = 0;
    J(3, 2)[8] = 1;
    J(3, 2)[9] = c40;
    J(3, 2)[10] = c41;
    J(3, 2)[11] = 0;
    J(3, 2)[12] = -c40;
    J(3, 2)[13] = -c41;
    J(3, 2)[14] = 0;
    J(3, 2)[15] = 0;
    J(3, 2)[16] = 0;
    J(3, 2)[17] = 1;
    J(3, 2)[18] = 0;
    J(3, 2)[19] = 0;
    J(3, 2)[20] = 0;
    J(3, 2)[21] = c41;
    J(3, 2)[22] = -c40;
    J(3, 2)[23] = 0;

    J(3, 3)[0] = -c41;
    J(3, 3)[1] = c40;
    J(3, 3)[2] = 0;
    J(3, 3)[3] = 0;
    J(3, 3)[4] = 0;
    J(3, 3)[5] = 0;
    J(3, 3)[6] = 0;
    J(3, 3)[7] = 0;
    J(3, 3)[8] = -1;
    J(3, 3)[9] = c40;
    J(3, 3)[10] = c41;
    J(3, 3)[11] = 0;
    J(3, 3)[12] = c40;
    J(3, 3)[13] = c41;
    J(3, 3)[14] = 0;
    J(3, 3)[15] = 0;
    J(3, 3)[16] = 0;
    J(3, 3)[17] = 1;
    J(3, 3)[18] = 0;
    J(3, 3)[19] = 0;
    J(3, 3)[20] = 0;
    J(3, 3)[21] = -c41;
    J(3, 3)[22] = c40;
    J(3, 3)[23] = 0;

  //math_policy::block6x4_type & WJTblock = WJT(0,0);
  //typename math_policy::real_type c2 = c1/20.0f;
    typename math_policy::real_type c3 = c1 / 20.0f;
    c30 /= 20.0f;
    c31 /= 20.0f;
    c40 /= 20.0f;
    c41 /= 20.0f;
    WJT(0, 0)[0] = 0;
    WJT(0, 0)[1] = 0;
    WJT(0, 0)[2] = 0;
    WJT(0, 0)[3] = 0;
    WJT(0, 0)[4] = 0;
    WJT(0, 0)[5] = 0;
    WJT(0, 0)[6] = 0;
    WJT(0, 0)[7] = 0;
    WJT(0, 0)[8] = 0;
    WJT(0, 0)[9] = 0;
    WJT(0, 0)[10] = 0;
    WJT(0, 0)[11] = 0;
    WJT(0, 0)[12] = 0;
    WJT(0, 0)[13] = 0;
    WJT(0, 0)[14] = 0;
    WJT(0, 0)[15] = 0;
    WJT(0, 0)[16] = 0;
    WJT(0, 0)[17] = 0;
    WJT(0, 0)[18] = 0;
    WJT(0, 0)[19] = 0;
    WJT(0, 0)[20] = 0;
    WJT(0, 0)[21] = 0;
    WJT(0, 0)[22] = 0;
    WJT(0, 0)[23] = 0;
  //WJTblock = WJT(1,0);
    WJT(1, 0)[0] = 0;
    WJT(1, 0)[1] = 0.05;
    WJT(1, 0)[2] = 0;
    WJT(1, 0)[3] = 0;
    WJT(1, 0)[4] = 0.05;
    WJT(1, 0)[5] = 0;
    WJT(1, 0)[6] = 0;
    WJT(1, 0)[7] = 0;
    WJT(1, 0)[8] = 0;
    WJT(1, 0)[9] = 0;
    WJT(1, 0)[10] = 0.05;
    WJT(1, 0)[11] = 0;
    WJT(1, 0)[12] = 0;
    WJT(1, 0)[13] = 0;
    WJT(1, 0)[14] = -0.05;
    WJT(1, 0)[15] = 0;
    WJT(1, 0)[16] = 0;
    WJT(1, 0)[17] = 0;
    WJT(1, 0)[18] = 0;
    WJT(1, 0)[19] = 0.05;
    WJT(1, 0)[20] = 0;
    WJT(1, 0)[21] = 0.05;
    WJT(1, 0)[22] = 0;
    WJT(1, 0)[23] = 0;
  //math_policy::block6x4_type & WJTblock = WJT(1,1);
    WJT(1, 1)[0] = -c3;
    WJT(1, 1)[1] = -c3;
    WJT(1, 1)[2] = 0;
    WJT(1, 1)[3] = 0;
    WJT(1, 1)[4] = -c3;
    WJT(1, 1)[5] = c3;
    WJT(1, 1)[6] = 0;
    WJT(1, 1)[7] = 0;
    WJT(1, 1)[8] = 0;
    WJT(1, 1)[9] = 0;
    WJT(1, 1)[10] = -0.05;
    WJT(1, 1)[11] = 0;
    WJT(1, 1)[12] = 0;
    WJT(1, 1)[13] = 0;
    WJT(1, 1)[14] = -c3;
    WJT(1, 1)[15] = -c3;
    WJT(1, 1)[16] = 0;
    WJT(1, 1)[17] = 0;
    WJT(1, 1)[18] = c3;
    WJT(1, 1)[19] = -c3;
    WJT(1, 1)[20] = 0;
    WJT(1, 1)[21] = 0.05;
    WJT(1, 1)[22] = 0;
    WJT(1, 1)[23] = 0;
  //math_policy::block6x4_type & WJTblock = WJT(1,2);
    WJT(1, 2)[0] = c30;
    WJT(1, 2)[1] = -c31;
    WJT(1, 2)[2] = 0;
    WJT(1, 2)[3] = 0;
    WJT(1, 2)[4] = -c31;
    WJT(1, 2)[5] = -c30;
    WJT(1, 2)[6] = 0;
    WJT(1, 2)[7] = 0;
    WJT(1, 2)[8] = 0;
    WJT(1, 2)[9] = 0;
    WJT(1, 2)[10] = -0.05;
    WJT(1, 2)[11] = 0;
    WJT(1, 2)[12] = 0;
    WJT(1, 2)[13] = 0;
    WJT(1, 2)[14] = -c31;
    WJT(1, 2)[15] = -c30;
    WJT(1, 2)[16] = 0;
    WJT(1, 2)[17] = 0;
    WJT(1, 2)[18] = -c30;
    WJT(1, 2)[19] = -c31;
    WJT(1, 2)[20] = 0;
    WJT(1, 2)[21] = 0.05;
    WJT(1, 2)[22] = 0;
    WJT(1, 2)[23] = 0;
  //WJTblock = WJT(2,1);
    WJT(2, 1)[0] = c3;
    WJT(2, 1)[1] = c3;
    WJT(2, 1)[2] = 0;
    WJT(2, 1)[3] = 0;
    WJT(2, 1)[4] = c3;
    WJT(2, 1)[5] = -c3;
    WJT(2, 1)[6] = 0;
    WJT(2, 1)[7] = 0;
    WJT(2, 1)[8] = 0;
    WJT(2, 1)[9] = 0;
    WJT(2, 1)[10] = 0.05;
    WJT(2, 1)[11] = 0;
    WJT(2, 1)[12] = 0;
    WJT(2, 1)[13] = 0;
    WJT(2, 1)[14] = c3;
    WJT(2, 1)[15] = c3;
    WJT(2, 1)[16] = 0;
    WJT(2, 1)[17] = 0;
    WJT(2, 1)[18] = c3;
    WJT(2, 1)[19] = c3;
    WJT(2, 1)[20] = 0;
    WJT(2, 1)[21] = 0.05;
    WJT(2, 1)[22] = 0;
    WJT(2, 1)[23] = 0;
  //math_policy::block6x4_type & WJTblock = WJT(2,3);
    WJT(2, 3)[0] = c31;
    WJT(2, 3)[1] = 0;
    WJT(2, 3)[2] = -c30;
    WJT(2, 3)[3] = 0;
    WJT(2, 3)[4] = -c30;
    WJT(2, 3)[5] = 0;
    WJT(2, 3)[6] = -c31;
    WJT(2, 3)[7] = 0;
    WJT(2, 3)[8] = 0;
    WJT(2, 3)[9] = 0.05;
    WJT(2, 3)[10] = 0;
    WJT(2, 3)[11] = 0;
    WJT(2, 3)[12] = 0;
    WJT(2, 3)[13] = c30;
    WJT(2, 3)[14] = 0;
    WJT(2, 3)[15] = c31;
    WJT(2, 3)[16] = 0;
    WJT(2, 3)[17] = c31;
    WJT(2, 3)[18] = 0;
    WJT(2, 3)[19] = -c30;
    WJT(2, 3)[20] = 0;
    WJT(2, 3)[21] = 0;
    WJT(2, 3)[22] = 0.05;
    WJT(2, 3)[23] = 0;

  //math_policy::block6x4_type & WJTblock = WJT(3,2);
    WJT(3, 2)[0] = -c30;
    WJT(3, 2)[1] = c31;
    WJT(3, 2)[2] = 0;
    WJT(3, 2)[3] = 0;
    WJT(3, 2)[4] = c31;
    WJT(3, 2)[5] = c30;
    WJT(3, 2)[6] = 0;
    WJT(3, 2)[7] = 0;
    WJT(3, 2)[8] = 0;
    WJT(3, 2)[9] = 0;
    WJT(3, 2)[10] = 0.05;
    WJT(3, 2)[11] = 0;
    WJT(3, 2)[12] = 0;
    WJT(3, 2)[13] = 0;
    WJT(3, 2)[14] = -c31;
    WJT(3, 2)[15] = -c30;
    WJT(3, 2)[16] = 0;
    WJT(3, 2)[17] = 0;
    WJT(3, 2)[18] = -c30;
    WJT(3, 2)[19] = c31;
    WJT(3, 2)[20] = 0;
    WJT(3, 2)[21] = 0.05;
    WJT(3, 2)[22] = 0;
    WJT(3, 2)[23] = 0;

  //math_policy::block6x4_type & WJTblock = WJT(3,3);
    WJT(3, 3)[0] = -c31;
    WJT(3, 3)[1] = 0;
    WJT(3, 3)[2] = c30;
    WJT(3, 3)[3] = 0;
    WJT(3, 3)[4] = c30;
    WJT(3, 3)[5] = 0;
    WJT(3, 3)[6] = c31;
    WJT(3, 3)[7] = 0;
    WJT(3, 3)[8] = 0;
    WJT(3, 3)[9] = -0.05;
    WJT(3, 3)[10] = 0;
    WJT(3, 3)[11] = 0;
    WJT(3, 3)[12] = 0;
    WJT(3, 3)[13] = c30;
    WJT(3, 3)[14] = 0;
    WJT(3, 3)[15] = -c31;
    WJT(3, 3)[16] = 0;
    WJT(3, 3)[17] = c31;
    WJT(3, 3)[18] = 0;
    WJT(3, 3)[19] = c30;
    WJT(3, 3)[20] = 0;
    WJT(3, 3)[21] = 0;
    WJT(3, 3)[22] = 0.05;
    WJT(3, 3)[23] = 0;

    BOOST_CHECK(false);
  //  prox_solver( ... fix me .... );

    typename math_policy::real_type tolerance = 1e-1; // that is 0.1 percent, an error on the fourth decimal
    BOOST_CHECK_CLOSE(lambda(0)[0], 5.327f, tolerance);
    BOOST_CHECK_CLOSE(lambda(0)[1], 0.3471f, tolerance);
    BOOST_CHECK_EQUAL(lambda(0)[2], 0);
    BOOST_CHECK_EQUAL(lambda(0)[3], 0);

    BOOST_CHECK_CLOSE(lambda(1)[0], 1.410f, tolerance);
    BOOST_CHECK_CLOSE(lambda(1)[1], -0.6658f, tolerance);
    BOOST_CHECK_EQUAL(lambda(1)[2], 0);
    BOOST_CHECK_EQUAL(lambda(1)[3], 0);

    BOOST_CHECK_CLOSE(lambda(2)[0], 1.887f, tolerance);
    BOOST_CHECK_CLOSE(lambda(2)[1], 0.2860f, tolerance);
    BOOST_CHECK_EQUAL(lambda(2)[2], 0);
    BOOST_CHECK_EQUAL(lambda(2)[3], 0);

    BOOST_CHECK_EQUAL(lambda(3)[2], 0);
    BOOST_CHECK_EQUAL(lambda(3)[3], 0);
    BOOST_CHECK_EQUAL(lambda(3)[2], 0);
    BOOST_CHECK_EQUAL(lambda(3)[3], 0);
}
*/
template <typename T>
void run_4_body_prox_solver_test_eigen(T i_am_a_dummy = T(0))
{

    // Create vectors with appropriate sizes
    Eigen::VectorX<T> q(28); // 4 bodies * 7
    Eigen::VectorX<T> Wdth(24); // 4 bodies * 6
    Eigen::VectorX<T> u(24); // 4 bodies * 6
    Eigen::VectorX<T> lambda(16); // 4 contacts * 4
    Eigen::VectorX<T> g(16); // 4 contacts * 4
    Eigen::VectorX<T> e(16); // 4 contacts * 4
    Eigen::VectorX<T> mu(16); // 4 contacts * 4

    // Create sparse matrices
    Eigen::SparseMatrix<T> J(
        16, 24); // 16 rows (4 contacts * 4), 24 cols (4 bodies * 6)
    Eigen::SparseMatrix<T> WJT(
        24, 16); // 24 rows (4 bodies * 6), 16 cols (4 contacts * 4)

    // Solver parameters
    prox::SolverParams<T> params;

    // Filling the vectors
    set_block7(q, 0, {0, 0, 0, 1.0, 0, 0, 0});
    set_block7(q, 1, {0, 2.0, 0, 1.0, 0, 0, 0});
    set_block7(q, 2, {1.4142, 3.4142, 0, 1.0, 0, 0, 0});
    set_block7(q, 3, {-0.51764, 3.9319, 0, 1.0, 0, 0, 0});

    set_block6(u, 0, {0, 0, 0, 0, 0, 0});
    set_block6(u, 1, {0, 0, 0, 0, 0, 0});
    set_block6(u, 2, {0, 0, 0, 0, 0, 0});
    set_block6(u, 3, {0, 0, 0, 0, 0, 0});

    set_block6(Wdth, 0, {0, 0, 0, 0, 0, 0});
    set_block6(Wdth, 1, {0, -0.0981, 0, 0, 0, 0});
    set_block6(Wdth, 2, {0, -0.0981, 0, 0, 0, 0});
    set_block6(Wdth, 3, {0, -0.0981, 0, 0, 0, 0});

    set_block4(g, 0, {0, 0, 0, 0});
    set_block4(g, 1, {0, 0, 0, 0});
    set_block4(g, 2, {0, 0, 0, 0});
    set_block4(g, 3, {0, 0, 0, 0});

    set_block4(e, 0, {1, 0, 0, 0});
    set_block4(e, 1, {1, 0, 0, 0});
    set_block4(e, 2, {1, 0, 0, 0});
    set_block4(e, 3, {1, 0, 0, 0});

    set_block4(mu, 0, {0, 0.5, 0.5, 0.5});
    set_block4(mu, 1, {0, 0.5, 0.5, 0.5});
    set_block4(mu, 2, {0, 0.5, 0.5, 0.5});
    set_block4(mu, 3, {0, 0.5, 0.5, 0.5});

    // Filling the matrices using triplets
    std::vector<Eigen::Triplet<T>> J_triplets;
    std::vector<Eigen::Triplet<T>> WJT_triplets;

    // Precompute constants
    T c1 = 1.0f / std::sqrt(2.0f);
    T c30 = 0.258849f;
    T c31 = 0.965927f;
    T c40 = 0.258849f;
    T c41 = 0.965927f;

    // J(0, 0) - contact 0, body 0
    J_triplets.emplace_back(0, 0, 0);
    J_triplets.emplace_back(0, 1, -1);
    J_triplets.emplace_back(0, 2, 0);
    J_triplets.emplace_back(0, 3, 0);
    J_triplets.emplace_back(0, 4, 0);
    J_triplets.emplace_back(0, 5, 0);
    J_triplets.emplace_back(1, 0, -1);
    J_triplets.emplace_back(1, 1, 0);
    J_triplets.emplace_back(1, 2, 0);
    J_triplets.emplace_back(1, 3, 0);
    J_triplets.emplace_back(1, 4, 0);
    J_triplets.emplace_back(1, 5, 1);
    J_triplets.emplace_back(2, 0, 0);
    J_triplets.emplace_back(2, 1, 0);
    J_triplets.emplace_back(2, 2, -1);
    J_triplets.emplace_back(2, 3, -1);
    J_triplets.emplace_back(2, 4, 0);
    J_triplets.emplace_back(2, 5, 0);
    J_triplets.emplace_back(3, 0, 0);
    J_triplets.emplace_back(3, 1, 0);
    J_triplets.emplace_back(3, 2, 0);
    J_triplets.emplace_back(3, 3, 0);
    J_triplets.emplace_back(3, 4, -1);
    J_triplets.emplace_back(3, 5, 0);

    // J(0, 1) - contact 0, body 1
    J_triplets.emplace_back(0, 6, 0);
    J_triplets.emplace_back(0, 7, 1);
    J_triplets.emplace_back(0, 8, 0);
    J_triplets.emplace_back(0, 9, 0);
    J_triplets.emplace_back(0, 10, 0);
    J_triplets.emplace_back(0, 11, 0);
    J_triplets.emplace_back(1, 6, 1);
    J_triplets.emplace_back(1, 7, 0);
    J_triplets.emplace_back(1, 8, 0);
    J_triplets.emplace_back(1, 9, 0);
    J_triplets.emplace_back(1, 10, 0);
    J_triplets.emplace_back(1, 11, 1);
    J_triplets.emplace_back(2, 6, 0);
    J_triplets.emplace_back(2, 7, 0);
    J_triplets.emplace_back(2, 8, 1);
    J_triplets.emplace_back(2, 9, -1);
    J_triplets.emplace_back(2, 10, 0);
    J_triplets.emplace_back(2, 11, 0);
    J_triplets.emplace_back(3, 6, 0);
    J_triplets.emplace_back(3, 7, 0);
    J_triplets.emplace_back(3, 8, 0);
    J_triplets.emplace_back(3, 9, 0);
    J_triplets.emplace_back(3, 10, 1);
    J_triplets.emplace_back(3, 11, 0);

    // J(1, 1) - contact 1, body 1
    J_triplets.emplace_back(4, 6, -c1);
    J_triplets.emplace_back(4, 7, -c1);
    J_triplets.emplace_back(4, 8, 0);
    J_triplets.emplace_back(4, 9, 0);
    J_triplets.emplace_back(4, 10, 0);
    J_triplets.emplace_back(4, 11, 0);
    J_triplets.emplace_back(5, 6, -c1);
    J_triplets.emplace_back(5, 7, c1);
    J_triplets.emplace_back(5, 8, 0);
    J_triplets.emplace_back(5, 9, 0);
    J_triplets.emplace_back(5, 10, 0);
    J_triplets.emplace_back(5, 11, 1);
    J_triplets.emplace_back(6, 6, 0);
    J_triplets.emplace_back(6, 7, 0);
    J_triplets.emplace_back(6, 8, -1);
    J_triplets.emplace_back(6, 9, -c1);
    J_triplets.emplace_back(6, 10, c1);
    J_triplets.emplace_back(6, 11, 0);
    J_triplets.emplace_back(7, 6, 0);
    J_triplets.emplace_back(7, 7, 0);
    J_triplets.emplace_back(7, 8, 0);
    J_triplets.emplace_back(7, 9, 0);
    J_triplets.emplace_back(7, 10, 0);
    J_triplets.emplace_back(7, 11, -c1);
    J_triplets.emplace_back(7, 6, -c1);

    // J(1, 2) - contact 1, body 2
    J_triplets.emplace_back(4, 12, c1);
    J_triplets.emplace_back(4, 13, c1);
    J_triplets.emplace_back(4, 14, 0);
    J_triplets.emplace_back(4, 15, 0);
    J_triplets.emplace_back(4, 16, 0);
    J_triplets.emplace_back(4, 17, 0);
    J_triplets.emplace_back(5, 12, c1);
    J_triplets.emplace_back(5, 13, -c1);
    J_triplets.emplace_back(5, 14, 0);
    J_triplets.emplace_back(5, 15, 0);
    J_triplets.emplace_back(5, 16, 0);
    J_triplets.emplace_back(5, 17, 1);
    J_triplets.emplace_back(6, 12, 0);
    J_triplets.emplace_back(6, 13, 0);
    J_triplets.emplace_back(6, 14, 1);
    J_triplets.emplace_back(6, 15, -c1);
    J_triplets.emplace_back(6, 16, c1);
    J_triplets.emplace_back(6, 17, 0);
    J_triplets.emplace_back(7, 12, 0);
    J_triplets.emplace_back(7, 13, 0);
    J_triplets.emplace_back(7, 14, 0);
    J_triplets.emplace_back(7, 15, 0);
    J_triplets.emplace_back(7, 16, 0);
    J_triplets.emplace_back(7, 17, c1);
    J_triplets.emplace_back(7, 12, c1);

    // J(2, 1) - contact 2, body 1
    J_triplets.emplace_back(8, 6, c30);
    J_triplets.emplace_back(8, 7, -c31);
    J_triplets.emplace_back(8, 8, 0);
    J_triplets.emplace_back(8, 9, 0);
    J_triplets.emplace_back(8, 10, 0);
    J_triplets.emplace_back(8, 11, 0);
    J_triplets.emplace_back(9, 6, -c31);
    J_triplets.emplace_back(9, 7, -c30);
    J_triplets.emplace_back(9, 8, 0);
    J_triplets.emplace_back(9, 9, 0);
    J_triplets.emplace_back(9, 10, 0);
    J_triplets.emplace_back(9, 11, 1);
    J_triplets.emplace_back(10, 6, 0);
    J_triplets.emplace_back(10, 7, 0);
    J_triplets.emplace_back(10, 8, -1);
    J_triplets.emplace_back(10, 9, -c31);
    J_triplets.emplace_back(10, 10, -c30);
    J_triplets.emplace_back(10, 11, 0);
    J_triplets.emplace_back(11, 6, 0);
    J_triplets.emplace_back(11, 7, 0);
    J_triplets.emplace_back(11, 8, 0);
    J_triplets.emplace_back(11, 9, 0);
    J_triplets.emplace_back(11, 10, 0);
    J_triplets.emplace_back(11, 11, c30);
    J_triplets.emplace_back(11, 6, -c31);

    // J(2, 3) - contact 2, body 3
    J_triplets.emplace_back(8, 18, -c30);
    J_triplets.emplace_back(8, 19, c31);
    J_triplets.emplace_back(8, 20, 0);
    J_triplets.emplace_back(8, 21, 0);
    J_triplets.emplace_back(8, 22, 0);
    J_triplets.emplace_back(8, 23, 0);
    J_triplets.emplace_back(9, 18, c31);
    J_triplets.emplace_back(9, 19, c30);
    J_triplets.emplace_back(9, 20, 0);
    J_triplets.emplace_back(9, 21, 0);
    J_triplets.emplace_back(9, 22, 0);
    J_triplets.emplace_back(9, 23, 1);
    J_triplets.emplace_back(10, 18, 0);
    J_triplets.emplace_back(10, 19, 0);
    J_triplets.emplace_back(10, 20, 1);
    J_triplets.emplace_back(10, 21, -c31);
    J_triplets.emplace_back(10, 22, -c30);
    J_triplets.emplace_back(10, 23, 0);
    J_triplets.emplace_back(11, 18, 0);
    J_triplets.emplace_back(11, 19, 0);
    J_triplets.emplace_back(11, 20, 0);
    J_triplets.emplace_back(11, 21, 0);
    J_triplets.emplace_back(11, 22, 0);
    J_triplets.emplace_back(11, 23, -c30);
    J_triplets.emplace_back(11, 18, c31);

    // J(3, 2) - contact 3, body 2
    J_triplets.emplace_back(12, 12, c41);
    J_triplets.emplace_back(12, 13, -c40);
    J_triplets.emplace_back(12, 14, 0);
    J_triplets.emplace_back(12, 15, 0);
    J_triplets.emplace_back(12, 16, 0);
    J_triplets.emplace_back(12, 17, 0);
    J_triplets.emplace_back(13, 12, 0);
    J_triplets.emplace_back(13, 13, 0);
    J_triplets.emplace_back(13, 14, 1);
    J_triplets.emplace_back(13, 15, c40);
    J_triplets.emplace_back(13, 16, c41);
    J_triplets.emplace_back(13, 17, 0);
    J_triplets.emplace_back(14, 12, -c40);
    J_triplets.emplace_back(14, 13, -c41);
    J_triplets.emplace_back(14, 14, 0);
    J_triplets.emplace_back(14, 15, 0);
    J_triplets.emplace_back(14, 16, 0);
    J_triplets.emplace_back(14, 17, 1);
    J_triplets.emplace_back(15, 12, 0);
    J_triplets.emplace_back(15, 13, 0);
    J_triplets.emplace_back(15, 14, 0);
    J_triplets.emplace_back(15, 15, 0);
    J_triplets.emplace_back(15, 16, 0);
    J_triplets.emplace_back(15, 17, c41);
    J_triplets.emplace_back(15, 12, -c40);

    // J(3, 3) - contact 3, body 3
    J_triplets.emplace_back(12, 18, -c41);
    J_triplets.emplace_back(12, 19, c40);
    J_triplets.emplace_back(12, 20, 0);
    J_triplets.emplace_back(12, 21, 0);
    J_triplets.emplace_back(12, 22, 0);
    J_triplets.emplace_back(12, 23, 0);
    J_triplets.emplace_back(13, 18, 0);
    J_triplets.emplace_back(13, 19, 0);
    J_triplets.emplace_back(13, 20, -1);
    J_triplets.emplace_back(13, 21, c40);
    J_triplets.emplace_back(13, 22, c41);
    J_triplets.emplace_back(13, 23, 0);
    J_triplets.emplace_back(14, 18, c40);
    J_triplets.emplace_back(14, 19, c41);
    J_triplets.emplace_back(14, 20, 0);
    J_triplets.emplace_back(14, 21, 0);
    J_triplets.emplace_back(14, 22, 0);
    J_triplets.emplace_back(14, 23, 1);
    J_triplets.emplace_back(15, 18, 0);
    J_triplets.emplace_back(15, 19, 0);
    J_triplets.emplace_back(15, 20, 0);
    J_triplets.emplace_back(15, 21, 0);
    J_triplets.emplace_back(15, 22, 0);
    J_triplets.emplace_back(15, 23, -c41);
    J_triplets.emplace_back(15, 18, c40);

    // WJT(0, 0) - body 0, contact 0
    for (int i = 0; i < 24; ++i) { WJT_triplets.emplace_back(i % 6, i / 6, 0); }

    // WJT(1, 0) - body 1, contact 0
    WJT_triplets.emplace_back(0, 0, 0);
    WJT_triplets.emplace_back(1, 0, 0.05);
    WJT_triplets.emplace_back(2, 0, 0);
    WJT_triplets.emplace_back(3, 0, 0);
    WJT_triplets.emplace_back(4, 0, 0.05);
    WJT_triplets.emplace_back(5, 0, 0);
    WJT_triplets.emplace_back(0, 1, 0);
    WJT_triplets.emplace_back(1, 1, 0);
    WJT_triplets.emplace_back(2, 1, 0);
    WJT_triplets.emplace_back(3, 1, 0);
    WJT_triplets.emplace_back(4, 1, 0);
    WJT_triplets.emplace_back(5, 1, 0.05);
    WJT_triplets.emplace_back(0, 2, 0);
    WJT_triplets.emplace_back(1, 2, 0);
    WJT_triplets.emplace_back(2, 2, 0);
    WJT_triplets.emplace_back(3, 2, -0.05);
    WJT_triplets.emplace_back(4, 2, 0);
    WJT_triplets.emplace_back(5, 2, 0);
    WJT_triplets.emplace_back(0, 3, 0);
    WJT_triplets.emplace_back(1, 3, 0.05);
    WJT_triplets.emplace_back(2, 3, 0);
    WJT_triplets.emplace_back(3, 3, 0.05);
    WJT_triplets.emplace_back(4, 3, 0);
    WJT_triplets.emplace_back(5, 3, 0);

    // Precompute scaled constants for WJT
    T c3 = c1 / 20.0f;
    T c30_scaled = c30 / 20.0f;
    T c31_scaled = c31 / 20.0f;
    T c40_scaled = c40 / 20.0f;
    T c41_scaled = c41 / 20.0f;

    // WJT(1, 1) - body 1, contact 1
    WJT_triplets.emplace_back(6, 4, -c3);
    WJT_triplets.emplace_back(7, 4, -c3);
    WJT_triplets.emplace_back(8, 4, 0);
    WJT_triplets.emplace_back(9, 4, 0);
    WJT_triplets.emplace_back(10, 4, -c3);
    WJT_triplets.emplace_back(11, 4, c3);
    WJT_triplets.emplace_back(6, 5, 0);
    WJT_triplets.emplace_back(7, 5, 0);
    WJT_triplets.emplace_back(8, 5, 0);
    WJT_triplets.emplace_back(9, 5, 0);
    WJT_triplets.emplace_back(10, 5, -0.05);
    WJT_triplets.emplace_back(11, 5, 0);
    WJT_triplets.emplace_back(6, 6, 0);
    WJT_triplets.emplace_back(7, 6, 0);
    WJT_triplets.emplace_back(8, 6, 0);
    WJT_triplets.emplace_back(9, 6, 0);
    WJT_triplets.emplace_back(10, 6, -c3);
    WJT_triplets.emplace_back(11, 6, -c3);
    WJT_triplets.emplace_back(6, 7, 0);
    WJT_triplets.emplace_back(7, 7, 0);
    WJT_triplets.emplace_back(8, 7, c3);
    WJT_triplets.emplace_back(9, 7, -c3);
    WJT_triplets.emplace_back(10, 7, 0);
    WJT_triplets.emplace_back(11, 7, 0.05);
    WJT_triplets.emplace_back(11, 7, 0);

    // WJT(1, 2) - body 1, contact 2
    WJT_triplets.emplace_back(6, 8, c30_scaled);
    WJT_triplets.emplace_back(7, 8, -c31_scaled);
    WJT_triplets.emplace_back(8, 8, 0);
    WJT_triplets.emplace_back(9, 8, 0);
    WJT_triplets.emplace_back(10, 8, -c31_scaled);
    WJT_triplets.emplace_back(11, 8, -c30_scaled);
    WJT_triplets.emplace_back(6, 9, 0);
    WJT_triplets.emplace_back(7, 9, 0);
    WJT_triplets.emplace_back(8, 9, 0);
    WJT_triplets.emplace_back(9, 9, 0);
    WJT_triplets.emplace_back(10, 9, -0.05);
    WJT_triplets.emplace_back(11, 9, 0);
    WJT_triplets.emplace_back(6, 10, 0);
    WJT_triplets.emplace_back(7, 10, 0);
    WJT_triplets.emplace_back(8, 10, 0);
    WJT_triplets.emplace_back(9, 10, 0);
    WJT_triplets.emplace_back(10, 10, -c31_scaled);
    WJT_triplets.emplace_back(11, 10, -c30_scaled);
    WJT_triplets.emplace_back(6, 11, 0);
    WJT_triplets.emplace_back(7, 11, 0);
    WJT_triplets.emplace_back(8, 11, -c30_scaled);
    WJT_triplets.emplace_back(9, 11, -c31_scaled);
    WJT_triplets.emplace_back(10, 11, 0);
    WJT_triplets.emplace_back(11, 11, 0.05);
    WJT_triplets.emplace_back(11, 11, 0);

    // WJT(2, 1) - body 2, contact 1
    WJT_triplets.emplace_back(12, 4, c3);
    WJT_triplets.emplace_back(13, 4, c3);
    WJT_triplets.emplace_back(14, 4, 0);
    WJT_triplets.emplace_back(15, 4, 0);
    WJT_triplets.emplace_back(16, 4, c3);
    WJT_triplets.emplace_back(17, 4, -c3);
    WJT_triplets.emplace_back(12, 5, 0);
    WJT_triplets.emplace_back(13, 5, 0);
    WJT_triplets.emplace_back(14, 5, 0);
    WJT_triplets.emplace_back(15, 5, 0);
    WJT_triplets.emplace_back(16, 5, 0.05);
    WJT_triplets.emplace_back(17, 5, 0);
    WJT_triplets.emplace_back(12, 6, 0);
    WJT_triplets.emplace_back(13, 6, 0);
    WJT_triplets.emplace_back(14, 6, 0);
    WJT_triplets.emplace_back(15, 6, 0);
    WJT_triplets.emplace_back(16, 6, c3);
    WJT_triplets.emplace_back(17, 6, c3);
    WJT_triplets.emplace_back(12, 7, 0);
    WJT_triplets.emplace_back(13, 7, 0);
    WJT_triplets.emplace_back(14, 7, c3);
    WJT_triplets.emplace_back(15, 7, c3);
    WJT_triplets.emplace_back(16, 7, 0);
    WJT_triplets.emplace_back(17, 7, 0.05);
    WJT_triplets.emplace_back(17, 7, 0);

    // WJT(2, 3) - body 2, contact 3
    WJT_triplets.emplace_back(12, 12, c31_scaled);
    WJT_triplets.emplace_back(13, 12, 0);
    WJT_triplets.emplace_back(14, 12, -c40_scaled);
    WJT_triplets.emplace_back(15, 12, 0);
    WJT_triplets.emplace_back(16, 12, -c40_scaled);
    WJT_triplets.emplace_back(17, 12, 0);
    WJT_triplets.emplace_back(12, 13, 0);
    WJT_triplets.emplace_back(13, 13, 0);
    WJT_triplets.emplace_back(14, 13, 0);
    WJT_triplets.emplace_back(15, 13, 0.05);
    WJT_triplets.emplace_back(16, 13, 0);
    WJT_triplets.emplace_back(17, 13, 0);
    WJT_triplets.emplace_back(12, 14, 0);
    WJT_triplets.emplace_back(13, 14, c30_scaled);
    WJT_triplets.emplace_back(14, 14, 0);
    WJT_triplets.emplace_back(15, 14, c31_scaled);
    WJT_triplets.emplace_back(16, 14, 0);
    WJT_triplets.emplace_back(17, 14, c31_scaled);
    WJT_triplets.emplace_back(12, 15, 0);
    WJT_triplets.emplace_back(13, 15, -c30_scaled);
    WJT_triplets.emplace_back(14, 15, 0);
    WJT_triplets.emplace_back(15, 15, 0);
    WJT_triplets.emplace_back(16, 15, 0);
    WJT_triplets.emplace_back(17, 15, 0.05);
    WJT_triplets.emplace_back(17, 15, 0);

    // WJT(3, 2) - body 3, contact 2
    WJT_triplets.emplace_back(18, 8, -c30_scaled);
    WJT_triplets.emplace_back(19, 8, c31_scaled);
    WJT_triplets.emplace_back(20, 8, 0);
    WJT_triplets.emplace_back(21, 8, 0);
    WJT_triplets.emplace_back(22, 8, c31_scaled);
    WJT_triplets.emplace_back(23, 8, c30_scaled);
    WJT_triplets.emplace_back(18, 9, 0);
    WJT_triplets.emplace_back(19, 9, 0);
    WJT_triplets.emplace_back(20, 9, 0);
    WJT_triplets.emplace_back(21, 9, 0);
    WJT_triplets.emplace_back(22, 9, 0.05);
    WJT_triplets.emplace_back(23, 9, 0);
    WJT_triplets.emplace_back(18, 10, 0);
    WJT_triplets.emplace_back(19, 10, 0);
    WJT_triplets.emplace_back(20, 10, 0);
    WJT_triplets.emplace_back(21, 10, 0);
    WJT_triplets.emplace_back(22, 10, -c31_scaled);
    WJT_triplets.emplace_back(23, 10, -c30_scaled);
    WJT_triplets.emplace_back(18, 11, 0);
    WJT_triplets.emplace_back(19, 11, 0);
    WJT_triplets.emplace_back(20, 11, -c30_scaled);
    WJT_triplets.emplace_back(21, 11, c31_scaled);
    WJT_triplets.emplace_back(22, 11, 0);
    WJT_triplets.emplace_back(23, 11, 0.05);
    WJT_triplets.emplace_back(23, 11, 0);

    // WJT(3, 3) - body 3, contact 3
    WJT_triplets.emplace_back(18, 12, -c41_scaled);
    WJT_triplets.emplace_back(19, 12, 0);
    WJT_triplets.emplace_back(20, 12, c40_scaled);
    WJT_triplets.emplace_back(21, 12, 0);
    WJT_triplets.emplace_back(22, 12, c40_scaled);
    WJT_triplets.emplace_back(23, 12, 0);
    WJT_triplets.emplace_back(18, 13, 0);
    WJT_triplets.emplace_back(19, 13, 0);
    WJT_triplets.emplace_back(20, 13, 0);
    WJT_triplets.emplace_back(21, 13, -0.05);
    WJT_triplets.emplace_back(22, 13, 0);
    WJT_triplets.emplace_back(23, 13, 0);
    WJT_triplets.emplace_back(18, 14, 0);
    WJT_triplets.emplace_back(19, 14, c30_scaled);
    WJT_triplets.emplace_back(20, 14, 0);
    WJT_triplets.emplace_back(21, 14, -c41_scaled);
    WJT_triplets.emplace_back(22, 14, 0);
    WJT_triplets.emplace_back(23, 14, c41_scaled);
    WJT_triplets.emplace_back(18, 15, 0);
    WJT_triplets.emplace_back(19, 15, c30_scaled);
    WJT_triplets.emplace_back(20, 15, 0);
    WJT_triplets.emplace_back(21, 15, c40_scaled);
    WJT_triplets.emplace_back(22, 15, 0);
    WJT_triplets.emplace_back(23, 15, 0.05);
    WJT_triplets.emplace_back(23, 15, 0);

    // Set the matrices from triplets
    J.setFromTriplets(J_triplets.begin(), J_triplets.end());
    WJT.setFromTriplets(WJT_triplets.begin(), WJT_triplets.end());

    BOOST_CHECK(false);
    // prox_solver(J, WJT, u, e, g, mu, lambda, params);

    T tolerance = 1e-1; // that is 0.1 percent, an error on the fourth decimal
    BOOST_CHECK_CLOSE(lambda(0), 5.327f, tolerance); // lambda(0)[0]
    BOOST_CHECK_CLOSE(lambda(1), 0.3471f, tolerance); // lambda(0)[1]
    BOOST_CHECK_EQUAL(lambda(2), 0); // lambda(0)[2]
    BOOST_CHECK_EQUAL(lambda(3), 0); // lambda(0)[3]

    BOOST_CHECK_CLOSE(lambda(4), 1.410f, tolerance); // lambda(1)[0]
    BOOST_CHECK_CLOSE(lambda(5), -0.6658f, tolerance); // lambda(1)[1]
    BOOST_CHECK_EQUAL(lambda(6), 0); // lambda(1)[2]
    BOOST_CHECK_EQUAL(lambda(7), 0); // lambda(1)[3]

    BOOST_CHECK_CLOSE(lambda(8), 1.887f, tolerance); // lambda(2)[0]
    BOOST_CHECK_CLOSE(lambda(9), 0.2860f, tolerance); // lambda(2)[1]
    BOOST_CHECK_EQUAL(lambda(10), 0); // lambda(2)[2]
    BOOST_CHECK_EQUAL(lambda(11), 0); // lambda(2)[3]

    BOOST_CHECK_EQUAL(lambda(12), 0); // lambda(3)[0]
    BOOST_CHECK_EQUAL(lambda(13), 0); // lambda(3)[1]
    BOOST_CHECK_EQUAL(lambda(14), 0); // lambda(3)[2]
    BOOST_CHECK_EQUAL(lambda(15), 0); // lambda(3)[3]
}
BOOST_AUTO_TEST_SUITE(prox_solver);

BOOST_AUTO_TEST_CASE(jacobi_2_body)
{
    using T = float;
    run_2_body_prox_solver_test_eigen<T>(/*prox::jacobi_solver<T>*/);
}

BOOST_AUTO_TEST_CASE(jacobi_3_body)
{
    using T = float;
    //run_3_body_prox_solver_test<T>(prox::jacobi_solver<T>);
    run_3_body_prox_solver_test_eigen<T>(/*prox::jacobi_solver<T>*/);
}

BOOST_AUTO_TEST_CASE(jacobi_4_body)
{
    using T = float;
    //    run_4_body_prox_solver_test<T>(prox::jacobi_solver<T>);
    run_4_body_prox_solver_test_eigen<T>();
}

BOOST_AUTO_TEST_CASE(gauss_seidel_2_body)
{
    using T = float;
    //    run_2_body_prox_solver_test_eigen<T>(prox::gauss_seidel_solver<T>);
    run_2_body_prox_solver_test_eigen<T>();
}

BOOST_AUTO_TEST_CASE(gauss_seidel_3_body)
{
    using T = float;
    //    run_3_body_prox_solver_test<T>(prox::gauss_seidel_solver<T>);
    run_3_body_prox_solver_test_eigen<T>(/*prox::gauss_seidel_solver<T>*/);
}

BOOST_AUTO_TEST_CASE(gauss_seidel_4_body)
{
    using T = float;
    //    run_4_body_prox_solver_test<T>(prox::gauss_seidel_solver<T>);
    run_4_body_prox_solver_test_eigen<T>();
}

BOOST_AUTO_TEST_SUITE_END();
