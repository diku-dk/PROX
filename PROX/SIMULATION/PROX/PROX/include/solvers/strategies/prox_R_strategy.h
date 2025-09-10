#ifndef PROX_R_STRATEGY_H
#define PROX_R_STRATEGY_H

#include "prox_enums.h"
#include "eigenhelperfunctions.h"

namespace prox
{

template <typename T>
void rstrategy_Eigen(strategy_type type, const Eigen::SparseMatrix<T>& J,
                     const Eigen::SparseMatrix<T>& WJT,
                     std::vector<Eigen::Matrix<T, 4, 4>>& R,
                     std::vector<Eigen::Matrix<T, 4, 4>>& nu)
{
    size_t M = J.rows() / 4; // Number of contact points
    assert(J.rows() % 4 == 0 && "J.rows() must be divisible by 4");

    switch (type)
    {
    case local_strategy:
        {
            R.resize(M);
            nu.resize(M);

            // Compute diagonal blocks of J * WJT
            for (size_t i = 0; i < M; ++i)
            {
                Eigen::Matrix<T, 4, 4> A_b
                    = J.middleRows(4 * i, 4) * WJT.middleCols(4 * i, 4);

                R[i].setZero();
                R[i](0, 0) = 1 / A_b(0, 0);
                R[i](1, 1) = 1 / A_b(1, 1);
                R[i](2, 2) = 1 / A_b(2, 2);
                R[i](3, 3) = 1 / A_b(3, 3);

                nu[i].setZero();
                nu[i](0, 0) = 0.9;
                nu[i](1, 1) = 0.9;
                nu[i](2, 2) = 0.9;
                nu[i](3, 3) = 0.9;
            }
            break;
        }
    case global_strategy:
        {
            R.resize(M);
            nu.resize(M);

            for (size_t i = 0; i < M; ++i)
            {
                R[i].setZero();
                R[i](0, 0) = 5;
                R[i](1, 1) = 5;
                R[i](2, 2) = 5;
                R[i](3, 3) = 5;

                nu[i].setZero();
                nu[i](0, 0) = 0.1;
                nu[i](1, 1) = 0.1;
                nu[i](2, 2) = 0.1;
                nu[i](3, 3) = 0.1;
            }
            break;
        }
    case blocked_strategy:
        {
            R.resize(M);
            nu.resize(M);

            // Compute diagonal blocks of J * WJT
            for (size_t i = 0; i < M; ++i)
            {
                Eigen::Matrix<T, 4, 4> A_b
                    = J.middleRows(4 * i, 4) * WJT.middleCols(4 * i, 4);

                // Extract 3x3 friction submatrix
                Eigen::Matrix<T, 3, 3> Aff;
                Aff(0, 0) = A_b(1, 1);
                Aff(0, 1) = A_b(1, 2);
                Aff(0, 2) = A_b(1, 3);
                Aff(1, 0) = A_b(2, 1);
                Aff(1, 1) = A_b(2, 2);
                Aff(1, 2) = A_b(2, 3);
                Aff(2, 0) = A_b(3, 1);
                Aff(2, 1) = A_b(3, 2);
                Aff(2, 2) = A_b(3, 3);

                // Invert the 3x3 matrix
                Aff = Aff.inverse();

                R[i].setZero();
                R[i](0, 0) = 1 / A_b(0, 0);
                R[i](1, 1) = Aff(0, 0);
                R[i](1, 2) = Aff(0, 1);
                R[i](1, 3) = Aff(0, 2);
                R[i](2, 1) = Aff(1, 0);
                R[i](2, 2) = Aff(1, 1);
                R[i](2, 3) = Aff(1, 2);
                R[i](3, 1) = Aff(2, 0);
                R[i](3, 2) = Aff(2, 1);
                R[i](3, 3) = Aff(2, 2);

                nu[i].setIdentity();
            }
            break;
        }
    default: throw std::logic_error("rstrategy(): unknown strategy type");
    }
}

} //namespace prox

// PROX_R_STRATEGY_H
#endif
