#ifndef PROX_SOLVER_H
#define PROX_SOLVER_H

#include <solvers/strategies/prox_R_strategy.h>
#include <solvers/sub/prox_normal_sub_solver.h>
#include <solvers/sub/prox_friction_sub_solver.h>

#include <solvers/prox_solver_params.h>
#include <solvers/prox_gauss_seidel_solver.h>
#include <solvers/prox_jacobi_solver.h>
#include <prox_math.h>

namespace prox
{

template <typename T>
void run_solver(const CRMatrix<4, 6, T>& J, const CRMatrix<6, 4, T>& WJT, const NCVec4<T>& b,
                const NCVec4<T>& mu, NCVec4<T>& lambda, const SolverParams<T>& params)
{
    switch (params.solver())
    {
    case solver_type::gauss_seidel: return gauss_seidel_solver(J, WJT, b, mu, lambda, params);
    case prox::solver_type::jacobi: return jacobi_solver(J, WJT, b, mu, lambda, params);
    default:                        assert(0);
    }
}
} //namespace prox

// PROX_SOLVER_H
#endif
