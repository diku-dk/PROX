#ifndef PROX_R_STRATEGY_H
#define PROX_R_STRATEGY_H

#include "prox_enums.h"
#include "prox_math.h"
#include "tiny_math_types.h"

namespace prox
{

template <typename T>
void rstrategy(strategy_type type, const CRMatrix<4, 6, T>& J,
               const CRMatrix<6, 4, T>& WJT, DiagonalMatrix<4, T>& R,
               DiagonalMatrix<4, T>& nu)
{
    switch (type)
    {
    case local_strategy:
        {
            // Get problem size
            size_t const M = J.nrows();
            nu.resize(M);
            R.resize(M);

            // From theory we want to have regularization
            // parameter selected as
            //
            // R = diag(r)
            //
            // where r_i = 1 / (J W J^T)_{ii}
            // In case of divergence then for the i'th index
            //
            //  r_i =   c_i r_i    if A_ii <= sum_{j \neq i} A_ij
            //              r_i    else
            //
            // where 0<c_i < 1 is some reduction paramter.
            //
            //
            // Apparently even for simple stacks all entries are diagonal
            // dominant and in such systems no reduction will happen. Thus
            // divergence is the only possibility. Therefor we are using a
            // global reduction.

            DiagonalMatrix<4, T> A;
            A.resize(J.nrows());
            sparse::diag_of_prod(J, WJT, A);

            // 2009-08-04 Kenny: we could save storage by reusing R to hold A, and the do inplace updating of R.

            for (size_t i = 0u; i < M; ++i)
            {
                auto& A_b = A(i);
                auto& R_b = R(i);
                auto& nu_b = nu(i);

                R_b(0, 0) = 1 / A_b(0, 0);
                R_b(1, 1) = 1 / A_b(1, 1);
                R_b(2, 2) = 1 / A_b(2, 2);
                R_b(3, 3) = 1 / A_b(3, 3);

                //assert off diagonals are indeed zero

                nu_b(0, 0) = 0.9f;
                nu_b(1, 1) = 0.9f;
                nu_b(2, 2) = 0.9f;
                nu_b(3, 3) = 0.9f;
            }
        }
        return;
    case global_strategy:
        {
            // Get problem size
            size_t const M = J.nrows();
            nu.resize(M);
            R.resize(M);

            // R-factor is some large positive scalar value. This should be 1/(|alpha_max| + |alpha_min| )
            // The reduction parameter is set to some user-specified positive value less than one.
            for (size_t i = 0u; i < M; ++i)
            {
                auto& R_b = R(i);
                auto& nu_b = nu(i);

                for (size_t j = 0u; j < 4; ++j)
                {
                    R_b(j, j) = 5;
                    nu_b(j, j) = 0.1f;
                }
            }
        }
        return;
    case blocked_strategy:
        {

            // Get problem size
            size_t const M = J.nrows();
            nu.resize(M);
            R.resize(M);

            // 2009-08-04 Kenny: Possible optimization; use R to store diagonal of A, then do
            // the inverse operations in place to get R. I actually think SPARSE supports
            // inversion of diagonal matrices so maybe one can just write R = sparse::inverse(R)
            // or something?

            DiagonalMatrix<4, T> A;
            A.resize(J.nrows());
            sparse::diag_of_prod(J, WJT, A);

            //
            // In the thesis of Studer he uses a block-diagonal matrix R, where
            // R_{ii} = A_{ii}^{-1}
            //
            for (size_t i = 0u; i < M; ++i)
            {
                auto& A_b = A(i);
                auto& R_b = R(i);
                auto& nu_b = nu(i);

                EigenMatrix3<T> Aff;
                Aff(0, 0) = A_b(1, 1);
                Aff(0, 1) = A_b(1, 2);
                Aff(0, 2) = A_b(1, 3);
                Aff(1, 0) = A_b(2, 1);
                Aff(1, 1) = A_b(2, 2);
                Aff(1, 2) = A_b(2, 3);
                Aff(2, 0) = A_b(3, 1);
                Aff(2, 1) = A_b(3, 2);
                Aff(2, 2) = A_b(3, 3);

                Aff = Aff.inverse();

                assert(fabs(A_b(0, 0)) > 0
                       || !"block_R_strategy(): divide by zero!");

                R_b(0, 0) = 1 / A_b(0, 0);
                R_b(1, 1) = Aff(0, 0);
                R_b(1, 2) = Aff(0, 1);
                R_b(1, 3) = Aff(0, 2);
                R_b(2, 1) = Aff(1, 0);
                R_b(2, 2) = Aff(1, 1);
                R_b(2, 3) = Aff(1, 2);
                R_b(3, 1) = Aff(2, 0);
                R_b(3, 2) = Aff(2, 1);
                R_b(3, 3) = Aff(2, 2);

                nu_b(0, 0) = 1;
                nu_b(1, 1) = 1;
                nu_b(2, 2) = 1;
                nu_b(3, 3) = 1;
            }
        }
        return;
    default: throw std::logic_error("About as broken as you'd expect.");
    }
}

} //namespace prox

// PROX_R_STRATEGY_H
#endif
