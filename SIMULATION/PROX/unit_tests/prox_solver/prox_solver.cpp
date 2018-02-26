#include <sparse.h>
#include <sparse_fill.h>
#include <prox_math_policy.h>
#include <solvers/prox_jacobi_solver.h>
#include <solvers/prox_gauss_seidel_solver.h>

#define BOOST_AUTO_TEST_MAIN
#include <boost/test/auto_unit_test.hpp>
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
template<typename math_policy, typename prox_solver_functor>
void run_2_body_prox_solver_test(prox_solver_functor const & prox_solver)
{
  
  typename math_policy::vector7_type q(2);
  typename math_policy::vector6_type Wdth(2), u(2);
  typename math_policy::vector4_type lambda(1), g(1), e(1), mu(1);
  typename math_policy::compressed4x6_type J(1,2,2);
  typename math_policy::compressed6x4_type WJT(2,1,2);
  
  prox::SolverParams<math_policy> params;
  
  //filling the vectors
  q(0)[0] = 0;  q(0)[1] = 0;   q(0)[2] = 0;  q(0)[3] = 0; q(0)[4] = 0;  q(0)[5] = 0;  q(0)[6] = -1.0;
  q(1)[0] = 0;  q(1)[1] = 2.0; q(1)[2] = 0;  q(1)[3] = 0; q(1)[4] = 0;  q(1)[5] = 0;  q(1)[6] = -1.0;
  
  u(0)[0] = 0;  u(0)[1] = 0;  u(0)[2] = 0;  u(0)[3] = 0; u(0)[4] = 0;  u(0)[5] = 0;
  u(1)[0] = 0;  u(1)[1] = 0;  u(1)[2] = 0;  u(1)[3] = 0; u(1)[4] = 0;  u(1)[5] = 0;
  
  Wdth(0)[0] = 0;  Wdth(0)[1] = 0;       Wdth(0)[2] = 0;  Wdth(0)[3] = 0; Wdth(0)[4] = 0;  Wdth(0)[5] = 0;
  Wdth(1)[0] = 0;  Wdth(1)[1] = -0.0981; Wdth(1)[2] = 0;  Wdth(1)[3] = 0; Wdth(1)[4] = 0;  Wdth(1)[5] = 0;
  
  g(0)[0]  = 0;  g(0)[1]  = 0;   g(0)[2]  = 0;   g(0)[3]  = 0;
  e(0)[0]  = 1;  e(0)[1]  = 0;   e(0)[2]  = 0;   e(0)[3]  = 0;
  mu(0)[0] = 0;  mu(0)[1] = 0.5; mu(0)[2] = 0.5; mu(0)[3] = 0.5;
  
  //filling the matrices
  //math_policy::block4x6_type & Jblock = J(0,0);
  J(0,0)[0]  = 0;  J(0,0)[1]  = -1; J(0,0)[2]  = 0;  J(0,0)[3]  = 0;  J(0,0)[4]  = 0;  J(0,0)[5]  = 0; 
  J(0,0)[6]  = -1; J(0,0)[7]  = 0;  J(0,0)[8]  = 0;  J(0,0)[9]  = 0;  J(0,0)[10] = 0;  J(0,0)[11] = 1; 
  J(0,0)[12] = 0;  J(0,0)[13] = 0;  J(0,0)[14] = -1; J(0,0)[15] = -1; J(0,0)[16] = 0;  J(0,0)[17] = 0; 
  J(0,0)[18] = 0;  J(0,0)[19] = 0;  J(0,0)[20] = 0;  J(0,0)[21] = 0;  J(0,0)[22] = -1; J(0,0)[23] = 0; 
  
  //math_policy::block4x6_type & Jblock2 = J(0,1);
  J(0,1)[0]  = 0;  J(0,1)[1]  = 1;  J(0,1)[2]  = 0;  J(0,1)[3]  = 0;  J(0,1)[4]  = 0;  J(0,1)[5]  = 0; 
  J(0,1)[6]  = 1;  J(0,1)[7]  = 0;  J(0,1)[8]  = 0;  J(0,1)[9]  = 0;  J(0,1)[10] = 0;  J(0,1)[11] = 1; 
  J(0,1)[12] = 0;  J(0,1)[13] = 0;  J(0,1)[14] = 1;  J(0,1)[15] = -1; J(0,1)[16] = 0;  J(0,1)[17] = 0; 
  J(0,1)[18] = 0;  J(0,1)[19] = 0;  J(0,1)[20] = 0;  J(0,1)[21] = 0;  J(0,1)[22] = 1;  J(0,1)[23] = 0; 
  
  //math_policy::block6x4_type & WJTblock = WJT(0,0);
  WJT(0,0)[0]  = 0; WJT(0,0)[1]  = 0; WJT(0,0)[2]  = 0; WJT(0,0)[3]  = 0;
  WJT(0,0)[4]  = 0; WJT(0,0)[5]  = 0; WJT(0,0)[6]  = 0; WJT(0,0)[7]  = 0;
  WJT(0,0)[8]  = 0; WJT(0,0)[9]  = 0; WJT(0,0)[10] = 0; WJT(0,0)[11] = 0;
  WJT(0,0)[12] = 0; WJT(0,0)[13] = 0; WJT(0,0)[14] = 0; WJT(0,0)[15] = 0;
  WJT(0,0)[16] = 0; WJT(0,0)[17] = 0; WJT(0,0)[18] = 0; WJT(0,0)[19] = 0;
  WJT(0,0)[20] = 0; WJT(0,0)[21] = 0; WJT(0,0)[22] = 0; WJT(0,0)[23] = 0;
  //WJTblock = WJT(1,0);
  WJT(1,0)[0]  = 0;    WJT(1,0)[1]  = 0.05; WJT(1,0)[2]  = 0;     WJT(1,0)[3]  = 0;
  WJT(1,0)[4]  = 0.05; WJT(1,0)[5]  = 0;    WJT(1,0)[6]  = 0;     WJT(1,0)[7]  = 0;
  WJT(1,0)[8]  = 0;    WJT(1,0)[9]  = 0;    WJT(1,0)[10] = 0.05;  WJT(1,0)[11] = 0;
  WJT(1,0)[12] = 0;    WJT(1,0)[13] = 0;    WJT(1,0)[14] = -0.05; WJT(1,0)[15] = 0;
  WJT(1,0)[16] = 0;    WJT(1,0)[17] = 0;    WJT(1,0)[18] = 0;     WJT(1,0)[19] = 0.05;
  WJT(1,0)[20] = 0;    WJT(1,0)[21] = 0.05; WJT(1,0)[22] = 0;     WJT(1,0)[23] = 0;
  
  
  
  BOOST_CHECK(false);
  //  prox_solver( ... fix me .... );
  typename math_policy::real_type tolerance = 1e-1; // that is 0.1 percent, an error on the fourth decimal
  BOOST_CHECK_CLOSE( lambda(0)[0], 1.962f, tolerance );
  BOOST_CHECK_EQUAL( lambda(0)[1], 0 );
  BOOST_CHECK_EQUAL( lambda(0)[2], 0 );
  BOOST_CHECK_EQUAL( lambda(0)[3], 0 );
  
  // We now apply a small velocity in the direction of the x axis on the sphere. 
  // This should result in friction. If the speed is 0.25, the impulse of the sphere
  // is P = m*v = 20*0.25 = 5. The bound on the friction impulse is 
  // mu*l_n = 0.5*1.962 = 0.981.
  lambda.clear();
  u(0)[0] = 0.25;
  BOOST_CHECK(false);
  //  prox_solver( ... fix me .... );
  
  BOOST_CHECK_CLOSE( lambda(0)[0], 1.962f, tolerance );
  BOOST_CHECK_CLOSE( lambda(0)[1], 0.981f, tolerance );
  BOOST_CHECK_EQUAL( lambda(0)[2], 0 );
  BOOST_CHECK_EQUAL( lambda(0)[3], 0 );
  
  // Applying a small velocity in the direction of the z axis on the sphere, should
  // result in friction. If the speed is 0.01, the impulse is P = m*v = 20*0.01 = 0.2. 
  // This changes the bounds on the friction to (0.980, 0.0390)
  lambda.clear();
  u(0)[2] = 0.01;
  
  BOOST_CHECK(false);
  //  prox_solver( ... fix me .... );
  
  BOOST_CHECK_CLOSE( lambda(0)[0], 1.962f, tolerance );
  BOOST_CHECK_CLOSE( lambda(0)[1], 0.980f, tolerance );
  BOOST_CHECK_CLOSE( lambda(0)[2], 0.0392f, tolerance );
  BOOST_CHECK_EQUAL( lambda(0)[3], 0 );
  
  // We now apply an angular velocity around the y axis on the sphere, to test Contensou friction
  // Setting w_y = 0.25, yields an upper bound of mu l_n = 0.981 (we assume the 2/3 factor to be
  // included in mu)
  u(0)[0] = 0; u(0)[2] = 0; u(0)[4] = 0.25;
  lambda.clear();

  BOOST_CHECK(false);
  //  prox_solver( ... fix me .... );
  
  BOOST_CHECK_CLOSE( lambda(0)[0], 1.962f, tolerance );
  BOOST_CHECK_EQUAL( lambda(0)[1], 0 );
  BOOST_CHECK_EQUAL( lambda(0)[2], 0 );
  BOOST_CHECK_CLOSE( lambda(0)[3], 0.981f, tolerance );
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
template<typename math_policy, typename prox_solver_functor>
void run_3_body_prox_solver_test(prox_solver_functor const & prox_solver )
{
  typename math_policy::vector7_type q(3);
  typename math_policy::vector6_type Wdth(3), u(3);
  typename math_policy::vector4_type lambda(2), g(2), e(2), mu(2);
  typename math_policy::compressed4x6_type J(2,3,4);
  typename math_policy::compressed6x4_type WJT(3,2,4);
  prox::SolverParams<math_policy> params;
  
  //filling the vectors
  q(0)[0] = 0;    q(0)[1] = 0;     q(0)[2] = 0;  q(0)[3] = 0; q(0)[4] = 0;  q(0)[5] = 0;  q(0)[6] = -1.0;
  q(1)[0] = 0;    q(1)[1] = 2.0;   q(1)[2] = 0;  q(1)[3] = 0; q(1)[4] = 0;  q(1)[5] = 0;  q(1)[6] = -1.0;
  q(2)[0] = 1.415;q(2)[1] = 3.415; q(2)[2] = 0;  q(2)[3] = 0; q(1)[4] = 0;  q(2)[5] = 0;  q(2)[6] = -1.0;
  
  u(0)[0] = 0;  u(0)[1] = 0;  u(0)[2] = 0;  u(0)[3] = 0; u(0)[4] = 0;  u(0)[5] = 0;
  u(1)[0] = 0;  u(1)[1] = 0;  u(1)[2] = 0;  u(1)[3] = 0; u(1)[4] = 0;  u(1)[5] = 0;
  u(2)[0] = 0;  u(2)[1] = 0;  u(2)[2] = 0;  u(2)[3] = 0; u(2)[4] = 0;  u(2)[5] = 0;
  
  Wdth(0)[0] = 0;  Wdth(0)[1] = 0;       Wdth(0)[2] = 0;  Wdth(0)[3] = 0; Wdth(0)[4] = 0;  Wdth(0)[5] = 0;
  Wdth(1)[0] = 0;  Wdth(1)[1] = -0.0981; Wdth(1)[2] = 0;  Wdth(1)[3] = 0; Wdth(1)[4] = 0;  Wdth(1)[5] = 0;
  Wdth(2)[0] = 0;  Wdth(2)[1] = -0.0981; Wdth(1)[2] = 0;  Wdth(2)[3] = 0; Wdth(2)[4] = 0;  Wdth(2)[5] = 0;
  
  g(0)[0]  = 0;  g(0)[1]  = 0;   g(0)[2]  = 0;   g(0)[3]  = 0;
  g(1)[0]  = 0;  g(1)[1]  = 0;   g(1)[2]  = 0;   g(1)[3]  = 0;
  e(0)[0]  = 1;  e(0)[1]  = 0;   e(0)[2]  = 0;   e(0)[3]  = 0;
  e(1)[0]  = 1;  e(1)[1]  = 0;   e(1)[2]  = 0;   e(1)[3]  = 0;
  mu(0)[0] = 0;  mu(0)[1] = 0.5; mu(0)[2] = 0.5; mu(0)[3] = 0.5;
  mu(1)[0] = 0;  mu(1)[1] = 0.5; mu(1)[2] = 0.5; mu(1)[3] = 0.5;
  
  //filling the matrices
  //math_policy::block4x6_type & Jblock = J(0,0);
  J(0,0)[0]  = 0;  J(0,0)[1]  = -1; J(0,0)[2]  = 0;  J(0,0)[3]  = 0;  J(0,0)[4]  = 0;  J(0,0)[5]  = 0; 
  J(0,0)[6]  = -1; J(0,0)[7]  = 0;  J(0,0)[8]  = 0;  J(0,0)[9]  = 0;  J(0,0)[10] = 0;  J(0,0)[11] = 1; 
  J(0,0)[12] = 0;  J(0,0)[13] = 0;  J(0,0)[14] = -1; J(0,0)[15] = -1; J(0,0)[16] = 0;  J(0,0)[17] = 0; 
  J(0,0)[18] = 0;  J(0,0)[19] = 0;  J(0,0)[20] = 0;  J(0,0)[21] = 0;  J(0,0)[22] = -1; J(0,0)[23] = 0; 
  
  //math_policy::block4x6_type & Jblock2 = J(0,1);
  J(0,1)[0]  = 0;  J(0,1)[1]  = 1;  J(0,1)[2]  = 0;  J(0,1)[3]  = 0;  J(0,1)[4]  = 0;  J(0,1)[5]  = 0; 
  J(0,1)[6]  = 1;  J(0,1)[7]  = 0;  J(0,1)[8]  = 0;  J(0,1)[9]  = 0;  J(0,1)[10] = 0;  J(0,1)[11] = 1; 
  J(0,1)[12] = 0;  J(0,1)[13] = 0;  J(0,1)[14] = 1;  J(0,1)[15] = -1; J(0,1)[16] = 0;  J(0,1)[17] = 0; 
  J(0,1)[18] = 0;  J(0,1)[19] = 0;  J(0,1)[20] = 0;  J(0,1)[21] = 0;  J(0,1)[22] = 1;  J(0,1)[23] = 0; 
  
  //math_policy::block4x6_type & Jblock = J(1,1);
  //math_policy::real_type c0 = sqrt(2.0f);
  typename math_policy::real_type c1 = 1.0f/sqrt(2.0f);
  J(1,1)[0]  = -c1; J(1,1)[1]  = -c1; J(1,1)[2]  = 0;  J(1,1)[3]  = 0;   J(1,1)[4]  = 0;   J(1,1)[5]  = 0; 
  J(1,1)[6]  = -c1; J(1,1)[7]  = c1;  J(1,1)[8]  = 0;  J(1,1)[9]  = 0;   J(1,1)[10] = 0;   J(1,1)[11] = 1; 
  J(1,1)[12] = 0;   J(1,1)[13] = 0;   J(1,1)[14] = -1; J(1,1)[15] = -c1; J(1,1)[16] = c1;  J(1,1)[17] = 0; 
  J(1,1)[18] = 0;   J(1,1)[19] = 0;   J(1,1)[20] = 0;  J(1,1)[21] = -c1; J(1,1)[22] = -c1; J(1,1)[23] = 0; 
  
  //math_policy::block4x6_type & Jblock2 = J(1,2);
  J(1,2)[0]  = c1;  J(1,2)[1]  = c1;  J(1,2)[2]  = 0;  J(1,2)[3]  = 0;   J(1,2)[4]  = 0;   J(1,2)[5]  = 0; 
  J(1,2)[6]  = c1;  J(1,2)[7]  = -c1; J(1,2)[8]  = 0;  J(1,2)[9]  = 0;   J(1,2)[10] = 0;   J(1,2)[11] = 1; 
  J(1,2)[12] = 0;   J(1,2)[13] = 0;   J(1,2)[14] = 1;  J(1,2)[15] = -c1; J(1,2)[16] = c1;  J(1,2)[17] = 0; 
  J(1,2)[18] = 0;   J(1,2)[19] = 0;   J(1,2)[20] = 0;  J(1,2)[21] = c1;  J(1,2)[22] = c1;  J(1,2)[23] = 0; 
  
  //math_policy::block6x4_type & WJTblock = WJT(0,0);
  WJT(0,0)[0]  = 0; WJT(0,0)[1]  = 0; WJT(0,0)[2]  = 0; WJT(0,0)[3]  = 0;
  WJT(0,0)[4]  = 0; WJT(0,0)[5]  = 0; WJT(0,0)[6]  = 0; WJT(0,0)[7]  = 0;
  WJT(0,0)[8]  = 0; WJT(0,0)[9]  = 0; WJT(0,0)[10] = 0; WJT(0,0)[11] = 0;
  WJT(0,0)[12] = 0; WJT(0,0)[13] = 0; WJT(0,0)[14] = 0; WJT(0,0)[15] = 0;
  WJT(0,0)[16] = 0; WJT(0,0)[17] = 0; WJT(0,0)[18] = 0; WJT(0,0)[19] = 0;
  WJT(0,0)[20] = 0; WJT(0,0)[21] = 0; WJT(0,0)[22] = 0; WJT(0,0)[23] = 0;
  //WJTblock = WJT(1,0);
  WJT(1,0)[0]  = 0;    WJT(1,0)[1]  = 0.05; WJT(1,0)[2]  = 0;     WJT(1,0)[3]  = 0;
  WJT(1,0)[4]  = 0.05; WJT(1,0)[5]  = 0;    WJT(1,0)[6]  = 0;     WJT(1,0)[7]  = 0;
  WJT(1,0)[8]  = 0;    WJT(1,0)[9]  = 0;    WJT(1,0)[10] = 0.05;  WJT(1,0)[11] = 0;
  WJT(1,0)[12] = 0;    WJT(1,0)[13] = 0;    WJT(1,0)[14] = -0.05; WJT(1,0)[15] = 0;
  WJT(1,0)[16] = 0;    WJT(1,0)[17] = 0;    WJT(1,0)[18] = 0;     WJT(1,0)[19] = 0.05;
  WJT(1,0)[20] = 0;    WJT(1,0)[21] = 0.05; WJT(1,0)[22] = 0;     WJT(1,0)[23] = 0;
  //math_policy::block6x4_type & WJTblock = WJT(1,1);
  //math_policy::real_type c2 = c0/20.0f;
  typename math_policy::real_type c3 = c1/20.0f;
  WJT(1,1)[0]  = -c3;  WJT(1,1)[1]  = -c3;   WJT(1,1)[2]  = 0;     WJT(1,1)[3]  = 0;
  WJT(1,1)[4]  = -c3;  WJT(1,1)[5]  = c3;    WJT(1,1)[6]  = 0;     WJT(1,1)[7]  = 0;
  WJT(1,1)[8]  = 0;    WJT(1,1)[9]  = 0;     WJT(1,1)[10] = -0.05; WJT(1,1)[11] = 0;
  WJT(1,1)[12] = 0;    WJT(1,1)[13] = 0;     WJT(1,1)[14] = -c3;   WJT(1,1)[15] = -c3;
  WJT(1,1)[16] = 0;    WJT(1,1)[17] = 0;     WJT(1,1)[18] = c3;    WJT(1,1)[19] = -c3;
  WJT(1,1)[20] = 0;    WJT(1,1)[21] = 0.05;  WJT(1,1)[22] = 0;     WJT(1,1)[23] = 0;
  //WJTblock = WJT(2,1);
  WJT(2,1)[0]  = c3;   WJT(2,1)[1]  = c3;    WJT(2,1)[2]  = 0;     WJT(2,1)[3]  = 0;
  WJT(2,1)[4]  = c3;   WJT(2,1)[5]  = -c3;   WJT(2,1)[6]  = 0;     WJT(2,1)[7]  = 0;
  WJT(2,1)[8]  = 0;    WJT(2,1)[9]  = 0;     WJT(2,1)[10] = 0.05;  WJT(2,1)[11] = 0;
  WJT(2,1)[12] = 0;    WJT(2,1)[13] = 0;     WJT(2,1)[14] = c3;    WJT(2,1)[15] = c3;
  WJT(2,1)[16] = 0;    WJT(2,1)[17] = 0;     WJT(2,1)[18] = c3;    WJT(2,1)[19] = c3;
  WJT(2,1)[20] = 0;    WJT(2,1)[21] = 0.05;  WJT(2,1)[22] = 0;     WJT(2,1)[23] = 0;
  
  BOOST_CHECK(false);
  //  prox_solver( ... fix me .... );

  typename math_policy::real_type tolerance = 1e-1; // that is 0.1 percent, an error on the fourth decimal
  BOOST_CHECK_CLOSE( lambda(0)[0], 3.428f, tolerance );
  BOOST_CHECK_CLOSE( lambda(0)[1], 0.5985f, tolerance );
  BOOST_CHECK_EQUAL( lambda(0)[2], 0 );
  BOOST_CHECK_EQUAL( lambda(0)[3], 0 );
  
  BOOST_CHECK_CLOSE( lambda(1)[0], 1.424f, tolerance );
  BOOST_CHECK_CLOSE( lambda(1)[1], -0.6498f, tolerance );
  BOOST_CHECK_EQUAL( lambda(1)[2], 0 );
  BOOST_CHECK_EQUAL( lambda(1)[3], 0 );
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
template<typename math_policy, typename prox_solver_functor>
void run_4_body_prox_solver_test(prox_solver_functor const & prox_solver)
{  
  typename math_policy::vector7_type q(4);
  typename math_policy::vector6_type Wdth(4), u(4);
  typename math_policy::vector4_type lambda(4), g(4), e(4), mu(4);
  typename math_policy::compressed4x6_type J(4,4,8);
  typename math_policy::compressed6x4_type WJT(4,4,8);
  prox::SolverParams<math_policy> params;
  
  //filling the vectors
  q(0)[0] = 0;      q(0)[1] = 0;      q(0)[2] = 0;  q(0)[3] = 1.0; q(0)[4] = 0;  q(0)[5] = 0;  q(0)[6] = 0;
  q(1)[0] = 0;      q(1)[1] = 2.0;    q(1)[2] = 0;  q(1)[3] = 1.0; q(1)[4] = 0;  q(1)[5] = 0;  q(1)[6] = 0;
  q(2)[0] = 1.4142; q(2)[1] = 3.4142; q(1)[2] = 0;  q(2)[3] = 1.0; q(2)[4] = 0;  q(2)[5] = 0;  q(2)[6] = 0;
  q(3)[0] = -0.51764; q(3)[1] = 3.9319;  q(3)[2] = 0;  q(3)[3] = 1.0; q(3)[4] = 0;  q(3)[5] = 0;  q(3)[6] = 0;
  
  u(0)[0] = 0;  u(0)[1] = 0;  u(0)[2] = 0;  u(0)[3] = 0; u(0)[4] = 0;  u(0)[5] = 0;
  u(1)[0] = 0;  u(1)[1] = 0;  u(1)[2] = 0;  u(1)[3] = 0; u(1)[4] = 0;  u(1)[5] = 0;
  u(2)[0] = 0;  u(2)[1] = 0;  u(2)[2] = 0;  u(2)[3] = 0; u(2)[4] = 0;  u(2)[5] = 0;
  u(3)[0] = 0;  u(3)[1] = 0;  u(3)[2] = 0;  u(3)[3] = 0; u(3)[4] = 0;  u(3)[5] = 0;
  
  Wdth(0)[0] = 0;  Wdth(0)[1] = 0;       Wdth(0)[2] = 0;  Wdth(0)[3] = 0; Wdth(0)[4] = 0;  Wdth(0)[5] = 0;
  Wdth(1)[0] = 0;  Wdth(1)[1] = -0.0981; Wdth(1)[2] = 0;  Wdth(1)[3] = 0; Wdth(1)[4] = 0;  Wdth(1)[5] = 0;
  Wdth(2)[0] = 0;  Wdth(2)[1] = -0.0981; Wdth(1)[2] = 0;  Wdth(2)[3] = 0; Wdth(2)[4] = 0;  Wdth(2)[5] = 0;
  Wdth(3)[0] = 0;  Wdth(3)[1] = -0.0981; Wdth(3)[2] = 0;  Wdth(3)[3] = 0; Wdth(3)[4] = 0;  Wdth(3)[5] = 0;
  
  g(0)[0]  = 0;  g(0)[1]  = 0;   g(0)[2]  = 0;   g(0)[3]  = 0;
  g(1)[0]  = 0;  g(1)[1]  = 0;   g(1)[2]  = 0;   g(1)[3]  = 0;
  g(2)[0]  = 0;  g(2)[1]  = 0;   g(2)[2]  = 0;   g(2)[3]  = 0;
  g(3)[0]  = 0;  g(3)[1]  = 0;   g(3)[2]  = 0;   g(3)[3]  = 0;
  
  e(0)[0]  = 1;  e(0)[1]  = 0;   e(0)[2]  = 0;   e(0)[3]  = 0;
  e(1)[0]  = 1;  e(1)[1]  = 0;   e(1)[2]  = 0;   e(1)[3]  = 0;
  e(2)[0]  = 1;  e(2)[1]  = 0;   e(2)[2]  = 0;   e(2)[3]  = 0;
  e(3)[0]  = 1;  e(3)[1]  = 0;   e(3)[2]  = 0;   e(3)[3]  = 0;
  
  mu(0)[0] = 0;  mu(0)[1] = 0.5; mu(0)[2] = 0.5; mu(0)[3] = 0.5;
  mu(1)[0] = 0;  mu(1)[1] = 0.5; mu(1)[2] = 0.5; mu(1)[3] = 0.5;
  mu(2)[0] = 0;  mu(2)[1] = 0.5; mu(2)[2] = 0.5; mu(2)[3] = 0.5;
  mu(3)[0] = 0;  mu(3)[1] = 0.5; mu(3)[2] = 0.5; mu(3)[3] = 0.5;
  
  //filling the matrices
  //math_policy::block4x6_type & Jblock = J(0,0);
  J(0,0)[0]  = 0;  J(0,0)[1]  = -1; J(0,0)[2]  = 0;  J(0,0)[3]  = 0;  J(0,0)[4]  = 0;  J(0,0)[5]  = 0; 
  J(0,0)[6]  = -1; J(0,0)[7]  = 0;  J(0,0)[8]  = 0;  J(0,0)[9]  = 0;  J(0,0)[10] = 0;  J(0,0)[11] = 1; 
  J(0,0)[12] = 0;  J(0,0)[13] = 0;  J(0,0)[14] = -1; J(0,0)[15] = -1; J(0,0)[16] = 0;  J(0,0)[17] = 0; 
  J(0,0)[18] = 0;  J(0,0)[19] = 0;  J(0,0)[20] = 0;  J(0,0)[21] = 0;  J(0,0)[22] = -1; J(0,0)[23] = 0; 
  
  //math_policy::block4x6_type & Jblock2 = J(0,1);
  J(0,1)[0]  = 0;  J(0,1)[1]  = 1;  J(0,1)[2]  = 0;  J(0,1)[3]  = 0;  J(0,1)[4]  = 0;  J(0,1)[5]  = 0; 
  J(0,1)[6]  = 1;  J(0,1)[7]  = 0;  J(0,1)[8]  = 0;  J(0,1)[9]  = 0;  J(0,1)[10] = 0;  J(0,1)[11] = 1; 
  J(0,1)[12] = 0;  J(0,1)[13] = 0;  J(0,1)[14] = 1;  J(0,1)[15] = -1; J(0,1)[16] = 0;  J(0,1)[17] = 0; 
  J(0,1)[18] = 0;  J(0,1)[19] = 0;  J(0,1)[20] = 0;  J(0,1)[21] = 0;  J(0,1)[22] = 1;  J(0,1)[23] = 0; 
  
  //math_policy::block4x6_type & Jblock = J(1,1);
  //math_policy::real_type c0 = sqrt(2.0f);
  typename math_policy::real_type c1 = 1.0f/sqrt(2.0f);
  J(1,1)[0]  = -c1; J(1,1)[1]  = -c1; J(1,1)[2]  = 0;  J(1,1)[3]  = 0;   J(1,1)[4]  = 0;   J(1,1)[5]  = 0; 
  J(1,1)[6]  = -c1; J(1,1)[7]  = c1;  J(1,1)[8]  = 0;  J(1,1)[9]  = 0;   J(1,1)[10] = 0;   J(1,1)[11] = 1; 
  J(1,1)[12] = 0;   J(1,1)[13] = 0;   J(1,1)[14] = -1; J(1,1)[15] = -c1; J(1,1)[16] = c1;  J(1,1)[17] = 0; 
  J(1,1)[18] = 0;   J(1,1)[19] = 0;   J(1,1)[20] = 0;  J(1,1)[21] = -c1;  J(1,1)[22] = -c1; J(1,1)[23] = 0; 
  
  //math_policy::block4x6_type & Jblock2 = J(1,2);
  J(1,2)[0]  = c1;  J(1,2)[1]  = c1;  J(1,2)[2]  = 0;  J(1,2)[3]  = 0;   J(1,2)[4]  = 0;   J(1,2)[5]  = 0; 
  J(1,2)[6]  = c1;  J(1,2)[7]  = -c1; J(1,2)[8]  = 0;  J(1,2)[9]  = 0;   J(1,2)[10] = 0;   J(1,2)[11] = 1; 
  J(1,2)[12] = 0;   J(1,2)[13] = 0;   J(1,2)[14] = 1;  J(1,2)[15] = -c1; J(1,2)[16] = c1;  J(1,2)[17] = 0; 
  J(1,2)[18] = 0;   J(1,2)[19] = 0;   J(1,2)[20] = 0;  J(1,2)[21] = c1;  J(1,2)[22] = c1;  J(1,2)[23] = 0; 
  
  // cp 3
  typename math_policy::real_type c30 = 0.258849;//0.259;
  typename math_policy::real_type c31 = 0.965927;//0.966;
  J(2,1)[0]  = c30;  J(2,1)[1]  = -c31; J(2,1)[2]  = 0;  J(2,1)[3]  = 0;    J(2,1)[4]  = 0;    J(2,1)[5]  = 0; 
  J(2,1)[6]  = -c31; J(2,1)[7]  = -c30; J(2,1)[8]  = 0;  J(2,1)[9]  = 0;    J(2,1)[10] = 0;    J(2,1)[11] = 1; 
  J(2,1)[12] = 0;    J(2,1)[13] = 0;    J(2,1)[14] = -1; J(2,1)[15] = -c31; J(2,1)[16] = -c30; J(2,1)[17] = 0; 
  J(2,1)[18] = 0;    J(2,1)[19] = 0;    J(2,1)[20] = 0;  J(2,1)[21] = c30;  J(2,1)[22] = -c31; J(2,1)[23] = 0; 
  
  //math_policy::block4x6_type & Jblock2 = J(1,2);
  J(2,3)[0]  = -c30; J(2,3)[1]  = c31;  J(2,3)[2]  = 0;  J(2,3)[3]  = 0;    J(2,3)[4]  = 0;    J(2,3)[5]  = 0; 
  J(2,3)[6]  = c31;  J(2,3)[7]  = c30;  J(2,3)[8]  = 0;  J(2,3)[9]  = 0;    J(2,3)[10] = 0;    J(2,3)[11] = 1; 
  J(2,3)[12] = 0;    J(2,3)[13] = 0;    J(2,3)[14] = 1;  J(2,3)[15] = -c31; J(2,3)[16] = -c30; J(2,3)[17] = 0; 
  J(2,3)[18] = 0;    J(2,3)[19] = 0;    J(2,3)[20] = 0;  J(2,3)[21] = -c30; J(2,3)[22] = c31;  J(2,3)[23] = 0; 
  
  // cp4
  typename math_policy::real_type c40 = 0.258849;//0.259;
  typename math_policy::real_type c41 = 0.965927;//0.966;
  J(3,2)[0]  = c41;  J(3,2)[1]  = -c40; J(3,2)[2]  = 0;  J(3,2)[3]  = 0;   J(3,2)[4]  = 0;    J(3,2)[5]  = 0; 
  J(3,2)[6]  = 0;    J(3,2)[7]  = 0;    J(3,2)[8]  = 1;  J(3,2)[9]  = c40; J(3,2)[10] = c41;  J(3,2)[11] = 0; 
  J(3,2)[12] = -c40; J(3,2)[13] = -c41; J(3,2)[14] = 0;  J(3,2)[15] = 0;   J(3,2)[16] = 0;    J(3,2)[17] = 1; 
  J(3,2)[18] = 0;    J(3,2)[19] = 0;    J(3,2)[20] = 0;  J(3,2)[21] = c41; J(3,2)[22] = -c40; J(3,2)[23] = 0; 
  
  J(3,3)[0]  = -c41; J(3,3)[1]  = c40;  J(3,3)[2]  = 0;  J(3,3)[3]  = 0;   J(3,3)[4]  = 0;    J(3,3)[5]  = 0; 
  J(3,3)[6]  = 0;    J(3,3)[7]  = 0;    J(3,3)[8]  = -1; J(3,3)[9]  = c40; J(3,3)[10] = c41;  J(3,3)[11] = 0; 
  J(3,3)[12] = c40;  J(3,3)[13] = c41;  J(3,3)[14] = 0;  J(3,3)[15] = 0;   J(3,3)[16] = 0;    J(3,3)[17] = 1; 
  J(3,3)[18] = 0;    J(3,3)[19] = 0;    J(3,3)[20] = 0;  J(3,3)[21] = -c41;J(3,3)[22] = c40;  J(3,3)[23] = 0; 
  
  //math_policy::block6x4_type & WJTblock = WJT(0,0);
  //typename math_policy::real_type c2 = c1/20.0f;
  typename math_policy::real_type c3 = c1/20.0f;
  c30 /= 20.0f;
  c31 /= 20.0f;
  c40 /= 20.0f;
  c41 /= 20.0f;
  WJT(0,0)[0]  = 0; WJT(0,0)[1]  = 0; WJT(0,0)[2]  = 0; WJT(0,0)[3]  = 0;
  WJT(0,0)[4]  = 0; WJT(0,0)[5]  = 0; WJT(0,0)[6]  = 0; WJT(0,0)[7]  = 0;
  WJT(0,0)[8]  = 0; WJT(0,0)[9]  = 0; WJT(0,0)[10] = 0; WJT(0,0)[11] = 0;
  WJT(0,0)[12] = 0; WJT(0,0)[13] = 0; WJT(0,0)[14] = 0; WJT(0,0)[15] = 0;
  WJT(0,0)[16] = 0; WJT(0,0)[17] = 0; WJT(0,0)[18] = 0; WJT(0,0)[19] = 0;
  WJT(0,0)[20] = 0; WJT(0,0)[21] = 0; WJT(0,0)[22] = 0; WJT(0,0)[23] = 0;
  //WJTblock = WJT(1,0);
  WJT(1,0)[0]  = 0;    WJT(1,0)[1]  = 0.05; WJT(1,0)[2]  = 0;     WJT(1,0)[3]  = 0;
  WJT(1,0)[4]  = 0.05; WJT(1,0)[5]  = 0;    WJT(1,0)[6]  = 0;     WJT(1,0)[7]  = 0;
  WJT(1,0)[8]  = 0;    WJT(1,0)[9]  = 0;    WJT(1,0)[10] = 0.05;  WJT(1,0)[11] = 0;
  WJT(1,0)[12] = 0;    WJT(1,0)[13] = 0;    WJT(1,0)[14] = -0.05; WJT(1,0)[15] = 0;
  WJT(1,0)[16] = 0;    WJT(1,0)[17] = 0;    WJT(1,0)[18] = 0;     WJT(1,0)[19] = 0.05;
  WJT(1,0)[20] = 0;    WJT(1,0)[21] = 0.05; WJT(1,0)[22] = 0;     WJT(1,0)[23] = 0;
  //math_policy::block6x4_type & WJTblock = WJT(1,1);
  WJT(1,1)[0]  = -c3;  WJT(1,1)[1]  = -c3;   WJT(1,1)[2]  = 0;     WJT(1,1)[3]  = 0;
  WJT(1,1)[4]  = -c3;  WJT(1,1)[5]  = c3;    WJT(1,1)[6]  = 0;     WJT(1,1)[7]  = 0;
  WJT(1,1)[8]  = 0;    WJT(1,1)[9]  = 0;     WJT(1,1)[10] = -0.05; WJT(1,1)[11] = 0;
  WJT(1,1)[12] = 0;    WJT(1,1)[13] = 0;     WJT(1,1)[14] = -c3;   WJT(1,1)[15] = -c3;
  WJT(1,1)[16] = 0;    WJT(1,1)[17] = 0;     WJT(1,1)[18] = c3;    WJT(1,1)[19] = -c3;
  WJT(1,1)[20] = 0;    WJT(1,1)[21] = 0.05;  WJT(1,1)[22] = 0;     WJT(1,1)[23] = 0;
  //math_policy::block6x4_type & WJTblock = WJT(1,2);
  WJT(1,2)[0]  = c30;  WJT(1,2)[1]  = -c31;  WJT(1,2)[2]  = 0;     WJT(1,2)[3]  = 0;
  WJT(1,2)[4]  = -c31; WJT(1,2)[5]  = -c30;  WJT(1,2)[6]  = 0;     WJT(1,2)[7]  = 0;
  WJT(1,2)[8]  = 0;    WJT(1,2)[9]  = 0;     WJT(1,2)[10] = -0.05; WJT(1,2)[11] = 0;
  WJT(1,2)[12] = 0;    WJT(1,2)[13] = 0;     WJT(1,2)[14] = -c31;  WJT(1,2)[15] = -c30;
  WJT(1,2)[16] = 0;    WJT(1,2)[17] = 0;     WJT(1,2)[18] = -c30;  WJT(1,2)[19] = -c31;
  WJT(1,2)[20] = 0;    WJT(1,2)[21] = 0.05;  WJT(1,2)[22] = 0;     WJT(1,2)[23] = 0;
  //WJTblock = WJT(2,1);
  WJT(2,1)[0]  = c3;   WJT(2,1)[1]  = c3;    WJT(2,1)[2]  = 0;     WJT(2,1)[3]  = 0;
  WJT(2,1)[4]  = c3;   WJT(2,1)[5]  = -c3;   WJT(2,1)[6]  = 0;     WJT(2,1)[7]  = 0;
  WJT(2,1)[8]  = 0;    WJT(2,1)[9]  = 0;     WJT(2,1)[10] = 0.05;  WJT(2,1)[11] = 0;
  WJT(2,1)[12] = 0;    WJT(2,1)[13] = 0;     WJT(2,1)[14] = c3;    WJT(2,1)[15] = c3;
  WJT(2,1)[16] = 0;    WJT(2,1)[17] = 0;     WJT(2,1)[18] = c3;    WJT(2,1)[19] = c3;
  WJT(2,1)[20] = 0;    WJT(2,1)[21] = 0.05;  WJT(2,1)[22] = 0;     WJT(2,1)[23] = 0;
  //math_policy::block6x4_type & WJTblock = WJT(2,3);
  WJT(2,3)[0]  = c31;  WJT(2,3)[1]  = 0;     WJT(2,3)[2]  = -c30;  WJT(2,3)[3]  = 0;
  WJT(2,3)[4]  = -c30; WJT(2,3)[5]  = 0;     WJT(2,3)[6]  = -c31;  WJT(2,3)[7]  = 0;
  WJT(2,3)[8]  = 0;    WJT(2,3)[9]  = 0.05;  WJT(2,3)[10] = 0;     WJT(2,3)[11] = 0;
  WJT(2,3)[12] = 0;    WJT(2,3)[13] = c30;   WJT(2,3)[14] = 0;     WJT(2,3)[15] = c31;
  WJT(2,3)[16] = 0;    WJT(2,3)[17] = c31;   WJT(2,3)[18] = 0;     WJT(2,3)[19] = -c30;
  WJT(2,3)[20] = 0;    WJT(2,3)[21] = 0;     WJT(2,3)[22] = 0.05;  WJT(2,3)[23] = 0;
  
  //math_policy::block6x4_type & WJTblock = WJT(3,2);
  WJT(3,2)[0]  = -c30; WJT(3,2)[1]  = c31;    WJT(3,2)[2]  = 0;    WJT(3,2)[3]  = 0;
  WJT(3,2)[4]  = c31;  WJT(3,2)[5]  = c30;    WJT(3,2)[6]  = 0;    WJT(3,2)[7]  = 0;
  WJT(3,2)[8]  = 0;    WJT(3,2)[9]  = 0;      WJT(3,2)[10] = 0.05; WJT(3,2)[11] = 0;
  WJT(3,2)[12] = 0;    WJT(3,2)[13] = 0;      WJT(3,2)[14] = -c31; WJT(3,2)[15] = -c30;
  WJT(3,2)[16] = 0;    WJT(3,2)[17] = 0;      WJT(3,2)[18] = -c30; WJT(3,2)[19] = c31;
  WJT(3,2)[20] = 0;    WJT(3,2)[21] = 0.05;   WJT(3,2)[22] = 0;    WJT(3,2)[23] = 0;
  
  //math_policy::block6x4_type & WJTblock = WJT(3,3);
  WJT(3,3)[0]  = -c31; WJT(3,3)[1]  = 0;     WJT(3,3)[2]  = c30;  WJT(3,3)[3]  = 0;
  WJT(3,3)[4]  = c30;  WJT(3,3)[5]  = 0;     WJT(3,3)[6]  = c31;  WJT(3,3)[7]  = 0;
  WJT(3,3)[8]  = 0;    WJT(3,3)[9]  = -0.05; WJT(3,3)[10] = 0;    WJT(3,3)[11] = 0;
  WJT(3,3)[12] = 0;    WJT(3,3)[13] = c30;   WJT(3,3)[14] = 0;    WJT(3,3)[15] = -c31;
  WJT(3,3)[16] = 0;    WJT(3,3)[17] = c31;   WJT(3,3)[18] = 0;    WJT(3,3)[19] = c30;
  WJT(3,3)[20] = 0;    WJT(3,3)[21] = 0;     WJT(3,3)[22] = 0.05; WJT(3,3)[23] = 0;
  
  BOOST_CHECK(false);
  //  prox_solver( ... fix me .... );
  
  typename math_policy::real_type tolerance = 1e-1; // that is 0.1 percent, an error on the fourth decimal
  BOOST_CHECK_CLOSE( lambda(0)[0], 5.327f, tolerance );
  BOOST_CHECK_CLOSE( lambda(0)[1], 0.3471f, tolerance );
  BOOST_CHECK_EQUAL( lambda(0)[2], 0 );
  BOOST_CHECK_EQUAL( lambda(0)[3], 0 );
  
  BOOST_CHECK_CLOSE( lambda(1)[0], 1.410f, tolerance );
  BOOST_CHECK_CLOSE( lambda(1)[1], -0.6658f, tolerance );
  BOOST_CHECK_EQUAL( lambda(1)[2], 0 );
  BOOST_CHECK_EQUAL( lambda(1)[3], 0 );
  
  BOOST_CHECK_CLOSE( lambda(2)[0], 1.887f, tolerance );
  BOOST_CHECK_CLOSE( lambda(2)[1], 0.2860f, tolerance );
  BOOST_CHECK_EQUAL( lambda(2)[2], 0 );
  BOOST_CHECK_EQUAL( lambda(2)[3], 0 );
  
  BOOST_CHECK_EQUAL( lambda(3)[2], 0 );
  BOOST_CHECK_EQUAL( lambda(3)[3], 0 );
  BOOST_CHECK_EQUAL( lambda(3)[2], 0 );
  BOOST_CHECK_EQUAL( lambda(3)[3], 0 );
}


BOOST_AUTO_TEST_SUITE(prox_solver);

BOOST_AUTO_TEST_CASE(jacobi_2_body)
{
  typedef prox::MathPolicy<float> math_policy;
  run_2_body_prox_solver_test<math_policy>( prox::jacobi_solver<math_policy> );
}

BOOST_AUTO_TEST_CASE(jacobi_3_body)
{
  typedef prox::MathPolicy<float> math_policy;
  run_3_body_prox_solver_test<math_policy>( prox::jacobi_solver<math_policy> );
}

BOOST_AUTO_TEST_CASE(jacobi_4_body)
{
  typedef prox::MathPolicy<float> math_policy;
  run_4_body_prox_solver_test<math_policy>( prox::jacobi_solver<math_policy> );
}

BOOST_AUTO_TEST_CASE(gauss_seidel_2_body)
{
  typedef prox::MathPolicy<float> math_policy;
  run_2_body_prox_solver_test<math_policy>( prox::gauss_seidel_solver<math_policy> );
}

BOOST_AUTO_TEST_CASE(gauss_seidel_3_body)
{
  typedef prox::MathPolicy<float> math_policy;
  run_3_body_prox_solver_test<math_policy>( prox::gauss_seidel_solver<math_policy> );
}

BOOST_AUTO_TEST_CASE(gauss_seidel_4_body)
{
  typedef prox::MathPolicy<float> math_policy;
  run_4_body_prox_solver_test<math_policy>( prox::gauss_seidel_solver<math_policy> );
}

BOOST_AUTO_TEST_SUITE_END();