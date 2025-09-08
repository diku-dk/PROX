#ifndef PROX_JACOBI_SOLVER_H
#define PROX_JACOBI_SOLVER_H

#include "prox_math.h"
#include <solvers/prox_solver_params.h>

#include <solvers/sub/prox_normal_sub_solver.h>
#include <solvers/sub/prox_friction_sub_solver.h>
#include <solvers/strategies/prox_R_strategy.h>

#include <util_profiling.h>
#include <util_log.h>

#include <sparse.h>

namespace prox
{

/**
   * This solver use the proximal map formulation of the constraints
   * on the motion of a rigid body system to compute the contact impulses
   * of the contact points. The system of constraints is solved using a
   * factorized Jacobi approach, where each constraint is solved atomically
   * in turn.
   */
template <typename T>
void jacobi_solver(const CRMatrix<4, 6, T>& J, const CRMatrix<6, 4, T>& WJT, const NCVec4<T>& b,
                   const NCVec4<T>& mu, NCVec4<T>& lambda, const SolverParams<T>& params)
{
    RECORD_VECTOR_NEW("convergence");
    RECORD_VECTOR_NEW("rfactor");

    START_TIMER("solver");

    size_t const K = J.nrows(); // Number of blocks

    size_t abs_conv_in_iteration
        = 0u; //used for profiling, to record in what iteration we found absolute convergence
    size_t rel_conv_in_iteration
        = 0u; //used for profiling, to record in what iteration we found relative convergence
    size_t count_divergence
        = 0u; //used for profiling, to record how many times we have discovered divergence

    //--- If warm starting is not used then clear the initial iterate
    if (!params.use_warm_starting()) { lambda.resize(K); }

    if (K == 0u) return;

    size_t in = 0; // Solution vector index that remembers the input iterate to the Jacobi scheme.
    size_t out = 1; // Solution vector index that remembers the output iterate to the Jacobi scheme.

    NCVec4<T> x
        [2]; // solution iterates, needed to do flip-flopping of solution vector in the Jacobi scheme
    x[in].resize(K);
    x[out].resize(K);

    if (params.use_warm_starting()) { x[in] = lambda; }

    NCVec4<T> residual;
    residual.resize(K);

    T last_residual_norm = std::numeric_limits<T>::max(); // Used to detect divergence.

    DiagonalMatrix<4, T> R, nu;

    rstrategy(params.r_factor_strategy(), J, WJT, R, nu);

    NCVec4<T> z;
    bool last_iteration_diverged = false;

    //--- Jacobi loops
    for (size_t iteration = 0u; iteration < params.max_iterations(); ++iteration)
    {
        RECORD_VECTOR_PUSH("rfactor", R(0)(1, 1));

        last_iteration_diverged = false;

        //--- Compute z = x - R(J W J^T x  + b)  = x - R( A x  + b)
        computeZ(x[in], R, J, WJT, b, z);

        for (size_t k = 0u; k < K; ++k)
        {
            const auto& mu_k = mu(k);
            const auto& z_k = z(k);
            const auto& x_k_in = x[in](k);
            auto& x_k_out = x[out](k);

            size_t const n = 0u;
            size_t const s = 1u;
            size_t const t = 2u;
            size_t const tau = 3u;

            normalSolver(params.normal_sub_solver(), z_k(n), x_k_out(n));

            frictionSolver(params.friction_sub_solver(), z_k(s), z_k(t),
                           z_k(tau), mu_k(s), mu_k(t), mu_k(tau), x_k_in(n),
                           x_k_out(s), x_k_out(t), x_k_out(tau));
        }

        //--- Compute residual, residual = lambda^k - lambda^(k+1)
        sparse::sub(x[in], x[out], residual);

        const auto residual_norm = computeInfNorm(residual);

        RECORD_VECTOR_PUSH("convergence", residual_norm);

        if (residual_norm < params.absolute_tolerance())
        {
            util::Log logging;

            logging << "jacobi_solver(): absolute convergence in " << iteration
                    << " iterations |residual| = " << residual_norm << util::Log::newline();

            abs_conv_in_iteration = iteration;

            break;
        }

        if (std::abs(residual_norm - last_residual_norm)
            < params.relative_tolerance() * last_residual_norm)
        {
            util::Log logging;

            logging << "jacobi_solver(): relative convergence in " << iteration << " iterations"
                    << util::Log::newline();

            rel_conv_in_iteration = iteration;

            break;
        }

        if (residual_norm > last_residual_norm)
        {
            util::Log logging;

            logging << "jacobi_solver(): divergence in " << iteration
                    << " iterations. |residual| = " << residual_norm << util::Log::newline();

            // Reduce R-factor and roll-back solution to last known good
            // iterate! (same as not doing a flip-flop on x-vectors).
            prod(nu, R);

            last_iteration_diverged = true;

            ++count_divergence;
        }
        else
        {
            last_residual_norm = residual_norm;
            std::swap(in, out);
            last_iteration_diverged = false;
        }
    }

    if (last_iteration_diverged)
        lambda = x[in];
    else
        lambda = x[out];

    RECORD("abs_conv", abs_conv_in_iteration);
    RECORD("rel_conv", rel_conv_in_iteration);
    RECORD("div_count", count_divergence);
    STOP_TIMER("solver");
}

template <typename T>
Eigen::Matrix<T, 4, 1> get_block4(const Eigen::VectorX<T>& v, size_t k)
{
    return Eigen::Matrix<T, 4, 1>(v(4 * k), v(4 * k + 1), v(4 * k + 2),
                                  v(4 * k + 3));
}

// Helper function to set a 4x1 block in a flat vector
template <typename T>
void set_block4(Eigen::VectorX<T>& v, size_t k,
                const Eigen::Matrix<T, 4, 1>& block)
{
    v(4 * k) = block(0);
    v(4 * k + 1) = block(1);
    v(4 * k + 2) = block(2);
    v(4 * k + 3) = block(3);
}

template <typename T>
void jacobi_solver_Eigen(const Eigen::SparseMatrix<T>& J,
                         const Eigen::SparseMatrix<T>& WJT,
                         const Eigen::VectorX<T>& b,
                         const Eigen::VectorX<T>& mu, Eigen::VectorX<T>& lambda,
                         const SolverParams<T>& params)
{
    RECORD_VECTOR_NEW("convergence");
    RECORD_VECTOR_NEW("rfactor");

    START_TIMER("solver");

    size_t const K_blocks = J.rows() / 4; // Number of contact points

    size_t abs_conv_in_iteration = 0;
    size_t rel_conv_in_iteration = 0;
    size_t count_divergence = 0;

    // If warm starting is not used then clear the initial iterate
    if (!params.use_warm_starting())
    {
        lambda.resize(4 * K_blocks);
        lambda.setZero();
    }

    if (K_blocks == 0) return;

    size_t in = 0; // Solution vector index for input iterate
    size_t out = 1; // Solution vector index for output iterate

    // Solution iterates for flip-flopping in Jacobi scheme
    Eigen::VectorX<T> x[2];
    x[in].resize(4 * K_blocks);
    x[out].resize(4 * K_blocks);

    if (params.use_warm_starting()) { x[in] = lambda; }
    else
    {
        x[in].setZero();
        x[out].setZero();
    }

    Eigen::VectorX<T> residual(4 * K_blocks);
    T last_residual_norm = std::numeric_limits<T>::max();

    std::vector<Eigen::Matrix<T, 4, 4>> R, nu;
    rstrategy(params.r_factor_strategy(), J, WJT, R, nu);

    Eigen::VectorX<T> z(4 * K_blocks);
    bool last_iteration_diverged = false;

    // Jacobi loops
    for (size_t iteration = 0; iteration < params.max_iterations(); ++iteration)
    {
        T val = (R.size() == 1) ? R[0](1, 1) : R[1](1, 1);
        RECORD_VECTOR_PUSH("rfactor", val);

        last_iteration_diverged = false;

        // Compute z = x - R(J W J^T x + b) = x - R(A x + b)
        // This would need a specialized computeZ function for flat vectors
        computeZ_flat(x[in], R, J, WJT, b, z);

        for (size_t k = 0; k < K_blocks; ++k)
        {
            Eigen::Matrix<T, 4, 1> mu_k = get_block4(mu, k);
            Eigen::Matrix<T, 4, 1> z_k = get_block4(z, k);
            Eigen::Matrix<T, 4, 1> x_k_in = get_block4(x[in], k);
            Eigen::Matrix<T, 4, 1> x_k_out = get_block4(x[out], k);

            size_t const n = 0;
            size_t const s = 1;
            size_t const t = 2;
            size_t const tau = 3;

            // Solve normal component
            normalSolver(params.normal_sub_solver(), z_k(n), x_k_out(n));

            // Solve friction components
            frictionSolver(params.friction_sub_solver(), z_k(s), z_k(t),
                           z_k(tau), mu_k(s), mu_k(t), mu_k(tau), x_k_in(n),
                           x_k_out(s), x_k_out(t), x_k_out(tau));

            // Store updated block back to flat vector
            set_block4(x[out], k, x_k_out);
        }

        // Compute residual: residual = x[in] - x[out]
        residual = x[in] - x[out];
        T residual_norm = computeInfNorm(residual);

        RECORD_VECTOR_PUSH("convergence", residual_norm);

        if (residual_norm < params.absolute_tolerance())
        {
            util::Log logging;
            logging << "jacobi_solver(): absolute convergence in " << iteration
                    << " iterations |residual| = " << residual_norm
                    << util::Log::newline();
            abs_conv_in_iteration = iteration;
            break;
        }

        if (std::abs(residual_norm - last_residual_norm)
            < params.relative_tolerance() * last_residual_norm)
        {
            util::Log logging;
            logging << "jacobi_solver(): relative convergence in " << iteration
                    << " iterations" << util::Log::newline();
            rel_conv_in_iteration = iteration;
            break;
        }

        if (residual_norm > last_residual_norm)
        {
            util::Log logging;
            logging << "jacobi_solver(): divergence in " << iteration
                    << " iterations. |residual| = " << residual_norm
                    << util::Log::newline();

            // Reduce R-factor and roll-back solution
            for (size_t i = 0; i < K_blocks; ++i) { R[i] = nu[i] * R[i]; }
            last_iteration_diverged = true;
            ++count_divergence;
        }
        else
        {
            last_residual_norm = residual_norm;
            std::swap(in, out);
            last_iteration_diverged = false;
        }
    }

    if (last_iteration_diverged) { lambda = x[in]; }
    else { lambda = x[out]; }

    RECORD("abs_conv", abs_conv_in_iteration);
    RECORD("rel_conv", rel_conv_in_iteration);
    RECORD("div_count", count_divergence);
    STOP_TIMER("solver");
}
} //namespace prox

// PROX_JACOBI_SOLVER_H
#endif
