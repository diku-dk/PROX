#ifndef PROX_GAUSS_SEIDEL_SOLVER_H
#define PROX_GAUSS_SEIDEL_SOLVER_H

#include "eigenhelperfunctions.h"
#include <solvers/sub/prox_normal_sub_solver.h>
#include <solvers/sub/prox_friction_sub_solver.h>
#include <solvers/strategies/prox_R_strategy.h>

#include <solvers/prox_solver_params.h>

#include <util_profiling.h>
#include <util_log.h>

namespace prox
{
/**
   * This solver use the proximal map formulation of the constraints
   * on the motion of a rigid body system to compute the contact impulses
   * of the contact points. The system of constraints is solved using a
   * factorized Gauss-Seidel approach, where each constraint is solved in
   * turn using the updated results from previously solved constraints.
   */

template <typename T>
void computeZk_Eigen(const Eigen::Matrix<T, 4, 1>& x_k,
                     const Eigen::VectorX<T>& w,
                     const Eigen::Matrix<T, 4, 4>& R_k,
                     const Eigen::SparseMatrix<T>& J,
                     const Eigen::Matrix<T, 4, 1>& b_k,
                     Eigen::Matrix<T, 4, 1>& z_k, size_t k)
{
    // Extract the k-th block row from J (rows 4*k to 4*k+3)
    Eigen::Matrix<T, 4, 1> Jw = J.middleRows(4 * k, 4) * w;

    // z_k = x_k - R_k * (Jw + b_k)
    z_k = x_k - R_k * (Jw + b_k);
}

template <typename T> T computeInfNorm_Eigen(const Eigen::VectorX<T>& x)
{
    assert(x.size() > 0 && "computeInfNorm(): zero sized vector");

    T max_val = 0;
    for (int i = 0; i < x.size(); ++i)
    {
        T abs_val = std::abs(x(i));
        if (abs_val > max_val) { max_val = abs_val; }
    }
    return max_val;
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
void gauss_seidel_solver_Eigen(const Eigen::SparseMatrix<T>& J,
                               const Eigen::SparseMatrix<T>& WJT,
                               const Eigen::VectorX<T>& b,
                               const Eigen::VectorX<T>& mu,
                               Eigen::VectorX<T>& lambda,
                               const SolverParams<T>& params)
{
    RECORD_VECTOR_NEW("convergence");
    RECORD_VECTOR_NEW("rfactor");

    START_TIMER("solver");

    size_t const K_blocks = J.rows() / 4; // Number of contact points

    size_t abs_conv_in_iteration = 0;
    size_t rel_conv_in_iteration = 0;
    size_t count_divergence = 0;

    if (!params.use_warm_starting())
    {
        lambda.resize(4 * K_blocks);
        lambda.setZero();
    }

    if (K_blocks == 0) return;

    Eigen::VectorX<T> x(4 * K_blocks); // Solution iterates

    if (params.use_warm_starting())
    {
        x = lambda; // Warm-starting
    }
    else { x.setZero(); }

    Eigen::VectorX<T> residual(4 * K_blocks);
    T last_residual_norm = std::numeric_limits<T>::max();

    std::vector<Eigen::Matrix<T, 4, 4>> R, nu;
    rstrategy_Eigen(params.r_factor_strategy(), J, WJT, R, nu);

    Eigen::VectorX<T> w(J.cols());
    w.setZero();

    T residual_norm;

    Eigen::SparseMatrix<T, Eigen::RowMajor> J_row = J; // convert once if needed
    J_row.makeCompressed(); // ensure compressed storage
    // WJT: keep as column-major (Eigen default). Make sure it's compressed:
    Eigen::SparseMatrix<T> WJT_col = WJT;
    WJT_col.makeCompressed();

    // Gauss-Seidel loops
    for (size_t iteration = 0; iteration < params.max_iterations(); ++iteration)
    {
        std::cerr << "Iteration " << iteration << "\n";
        // Safety guard for R factor
        T val = (R.size() == 1) ? R[0](1, 1) : R[1](1, 1);
        RECORD_VECTOR_PUSH("rfactor", val);

        // Loop over contact points
        for (size_t k = 0; k < K_blocks; ++k)
        {
            //std::cerr << "k: " << k << "/" << K_blocks << "\n";
            using Vec4 = Eigen::Matrix<T, 4, 1>;

            // map the 4-element blocks directly into the big vectors (no copy)
            Eigen::Map<Vec4> xk_map(x.data()
                                    + 4 * k); // references x[4*k..4*k+3]
            Eigen::Map<const Vec4> b_k_map(b.data() + 4 * k);
            Eigen::Map<const Vec4> mu_k_map(mu.data() + 4 * k);

            // save old x_k
            Vec4 old_xk = xk_map;

            //compute Jw = J_k * w  (J_k = rows 4*k .. 4*k+3)
            Vec4 Jw;
            Jw.setZero();
            for (int r = 0; r < 4; ++r)
            {
                int row = static_cast<int>(4 * k + r);
                for (typename Eigen::SparseMatrix<
                         T, Eigen::RowMajor>::InnerIterator it(J_row, row);
                     it; ++it)
                {
                    // it.col() is the column index; w[it.col()] is dense-vector access
                    Jw(r) += it.value() * w[it.col()];
                }
            }

            // z_k = x_k - R[k] * (Jw + b_k)
            Vec4 z_k = xk_map - R[k] * (Jw + b_k_map);

            // solve local (these are scalar calls that update xk_map in place)
            // note: normalSolver and frictionSolver expect scalars; adapt to how they modify x_k
            // Copy z_k components into temporaries if your solvers require references
            //--- Solve lambda_n = prox_{R^+}( lambda_n - r (A lambda_n + b))
            normalSolver(params.normal_sub_solver(), z_k(0), xk_map(0));
            //--- Solve lambda_f = prox_C( lambda_f - r (A lambda_f + b))
            frictionSolver(params.friction_sub_solver(), z_k(1), z_k(2), z_k(3),
                           mu_k_map(1), mu_k_map(2), mu_k_map(3), xk_map(0),
                           xk_map(1), xk_map(2), xk_map(3));

            // delta_x = x_k_new - old_xk
            Vec4 delta_x = xk_map - old_xk;

            // update w: w += WJT.blockColumns(4*k..4*k+3) * delta_x
            // efficient column-wise accumulation (WJT_col is column-major)
            int baseCol = 4 * static_cast<int>(k);
            for (int local_col = 0; local_col < 4; ++local_col)
            {
                int col = baseCol + local_col;
                T coeff = delta_x(local_col); // scalar
                if (coeff == T(0)) continue; // small cheap optimization
                for (typename Eigen::SparseMatrix<T>::InnerIterator it(WJT_col,
                                                                       col);
                     it; ++it)
                {
                    // it.row() gives row index into w
                    w[it.row()] += it.value() * coeff;
                }
            }
        }

        // Compute residual
        residual = lambda - x;
        residual_norm = computeInfNorm_Eigen(residual);

        RECORD_VECTOR_PUSH("convergence", residual_norm);

        if (residual_norm < params.absolute_tolerance())
        {
            util::Log logging;
            logging << "gauss_seidel_solver(): absolute convergence in "
                    << iteration << " iterations |residual| = " << residual_norm
                    << util::Log::newline();
            abs_conv_in_iteration = iteration;
            break;
        }

        if (std::abs(residual_norm - last_residual_norm)
            < params.relative_tolerance() * last_residual_norm)
        {
            util::Log logging;
            logging << "gauss_seidel_solver(): relative convergence in "
                    << iteration << " iterations" << util::Log::newline();
            rel_conv_in_iteration = iteration;
            break;
        }

        if (residual_norm > last_residual_norm)
        {
            util::Log logging;
            logging << "gauss_seidel_solver(): divergence in " << iteration
                    << " iterations. |residual| = " << residual_norm
                    << util::Log::newline();

            // Reduce R-factor and roll-back solution
            for (size_t i = 0; i < K_blocks; ++i) { R[i] = nu[i] * R[i]; }
            x = lambda;
            ++count_divergence;
        }
        else
        {
            last_residual_norm = residual_norm;
            lambda = x;
        }
    }

    lambda = x;

    RECORD("abs_conv", abs_conv_in_iteration);
    RECORD("rel_conv", rel_conv_in_iteration);
    RECORD("div_count", count_divergence);
    STOP_TIMER("solver");
}

} //namespace prox

// PROX_GAUSS_SEIDEL_SOLVER_H
#endif
