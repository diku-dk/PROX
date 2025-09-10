#ifndef PROX_GET_INVERSE_MASS_MATRIX_H
#define PROX_GET_INVERSE_MASS_MATRIX_H

#include "prox_math.h"
#include <prox_update_inertia_tensor.h>

#include <eigenhelperfunctions.h>

#include <cassert>

namespace prox
{

  // 2009-08-13 Kenny code reivew: Optimization replace diagonal6x6_type with diagonal_mass_type, maybe wait to optimize until all it working

/*template <typename T, typename Iterator>
inline void get_inverse_mass_matrix(Iterator begin, Iterator end, DiagonalMatrix<6, T>& W)
{
    size_t const N = std::distance(begin,end);

    W.resize( N );

    size_t index = 0u;
    for (auto body = begin; body != end; ++body, ++index)
    {
        T inv_mass = 0;
        EigenMatrix3<T> inv_I{
            {0, 0, 0},
            {0, 0, 0},
            {0, 0, 0}
        };

        if (!body->is_fixed() && !body->is_scripted())
        {
            assert(fabs(body->get_mass()) > 0 || !"get_inverse_mass_matrix(): Divide by zero!");

            inv_mass = 1 / body->get_mass();

            assert(is_number(inv_mass) || !"get_inverse_mass_matrix(): Nan");
            assert(is_finite(inv_mass) || !"get_inverse_mass_matrix(): Inf");
            assert(inv_mass > 0 || !"get_inverse_mass_matrix(): Negative mass");

            auto I_bf = body->get_inertia_bf();
            auto R = EigenMatrix3<T>(body->get_orientation());

            detail::update_inertia_tensor(R, I_bf, inv_I);

            inv_I = inv_I.inverse().eval();
      }

      auto& b = W(index);

      b(0,0) = inv_mass;
      b(1,1) = inv_mass;
      b(2,2) = inv_mass;
      b(3,3) = inv_I(0,0);
      b(3,4) = inv_I(0,1);
      b(3,5) = inv_I(0,2);
      b(4,3) = inv_I(1,0);
      b(4,4) = inv_I(1,1);
      b(4,5) = inv_I(1,2);
      b(5,3) = inv_I(2,0);
      b(5,4) = inv_I(2,1);
      b(5,5) = inv_I(2,2);
    }
}*/

template <typename T, typename Iterator>
inline void get_inverse_mass_matrix_eigen(Iterator begin, Iterator end,
                                          Eigen::SparseMatrix<T>& W)
{
    size_t const N = std::distance(begin, end);
    std::vector<Eigen::Triplet<T>> triplets;
    triplets.reserve(N * 36);

    size_t index = 0;
    for (auto body = begin; body != end; ++body, ++index)
    {
        T inv_mass = 0;
        Eigen::Matrix<T, 3, 3> inv_I = Eigen::Matrix<T, 3, 3>::Zero();

        if (!body->is_fixed() && !body->is_scripted())
        {
            assert(std::abs(body->get_mass()) > 0
                   && "get_inverse_mass_matrix(): Divide by zero!");

            inv_mass = 1 / body->get_mass();

            assert(is_number(inv_mass) && "get_inverse_mass_matrix(): Nan");
            assert(is_finite(inv_mass) && "get_inverse_mass_matrix(): Inf");
            assert(inv_mass > 0 && "get_inverse_mass_matrix(): Negative mass");

            // Replicate the original calculation exactly
            auto I_bf = body->get_inertia_bf();
            auto R = Eigen::Matrix<T, 3, 3>(
                body->get_orientation()); // Same as original

            // Use the same update_inertia_tensor function
            Eigen::Matrix<T, 3, 3> temp_I;
            detail::update_inertia_tensor(R, I_bf, temp_I);
            inv_I = temp_I.inverse();
        }

        size_t start_idx = index * 6;

        // Set the diagonal entries for mass
        triplets.emplace_back(start_idx + 0, start_idx + 0, inv_mass);
        triplets.emplace_back(start_idx + 1, start_idx + 1, inv_mass);
        triplets.emplace_back(start_idx + 2, start_idx + 2, inv_mass);

        // Set the entries for inertia tensor - ensure same order as original
        triplets.emplace_back(start_idx + 3, start_idx + 3, inv_I(0, 0));
        triplets.emplace_back(start_idx + 3, start_idx + 4, inv_I(0, 1));
        triplets.emplace_back(start_idx + 3, start_idx + 5, inv_I(0, 2));

        triplets.emplace_back(start_idx + 4, start_idx + 3, inv_I(1, 0));
        triplets.emplace_back(start_idx + 4, start_idx + 4, inv_I(1, 1));
        triplets.emplace_back(start_idx + 4, start_idx + 5, inv_I(1, 2));

        triplets.emplace_back(start_idx + 5, start_idx + 3, inv_I(2, 0));
        triplets.emplace_back(start_idx + 5, start_idx + 4, inv_I(2, 1));
        triplets.emplace_back(start_idx + 5, start_idx + 5, inv_I(2, 2));
    }

    W.resize(6 * N, 6 * N);
    W.setFromTriplets(triplets.begin(), triplets.end());
}

template <typename T>
void verify_inverse_mass_matrix(const DiagonalMatrix<6, T>& W_orig,
                                const Eigen::SparseMatrix<T>& W_eigen,
                                double tolerance = 1e-15)
{
    size_t N = W_orig.size();

    for (size_t i = 0; i < N; ++i)
    {
        const auto& block_orig = W_orig(i);
        size_t start_idx = i * 6;

        for (int row = 0; row < 6; ++row)
        {
            for (int col = 0; col < 6; ++col)
            {
                T val_orig = block_orig(row, col);
                T val_eigen = W_eigen.coeff(start_idx + row, start_idx + col);

                if (std::abs(val_orig - val_eigen) > tolerance)
                {
                    std::cout << "Mismatch at body " << i << " (" << row << ", "
                              << col << "): " << val_orig << " vs " << val_eigen
                              << std::endl;
                }
            }
        }
    }
}
} // namespace prox
// PROX_GET_INVERSE_MASS_MATRIX_H
#endif
