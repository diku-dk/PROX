#ifndef PROX_SOLVER_H
#define PROX_SOLVER_H

#include <solvers/strategies/prox_R_strategy.h>
#include <solvers/sub/prox_normal_sub_solver.h>
#include <solvers/sub/prox_friction_sub_solver.h>

#include <solvers/prox_solver_params.h>

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

/*
template <typename T>
struct SimpleSolver
{
    SimpleSolver(        logging << "bind_solver(): using jacobi solver"<< util::Log::newline();
                 return SolverBinder<M>( &jacobi_solver<M> );
                 
                 case gauss_seidel:
                 logging << "bind_solver(): using gauss seidel solver"<< util::Log::newline();
                 return SolverBinder<M>( &gauss_seidel_solver<M> );
)
*/
  /**
   * A solver functor.
   */
template <typename M> class NoSuchSolver
{
public:
    virtual void operator()(typename M::compressed4x6_type const&,
                            typename M::compressed6x4_type const&, typename M::vector4_type const&,
                            typename M::vector4_type const&, typename M::vector4_type&,
                            RStrategy<M> const&, NormalSubSolver<typename M::real_type> const&,
                            FrictionSubSolver<typename M::real_type> const&, SolverParams<M> const&,
                            M const&) const
        = 0;
};

} //namespace prox

// PROX_SOLVER_H
#endif
