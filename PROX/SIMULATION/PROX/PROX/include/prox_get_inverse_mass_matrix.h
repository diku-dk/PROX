#ifndef PROX_GET_INVERSE_MASS_MATRIX_H
#define PROX_GET_INVERSE_MASS_MATRIX_H

#include "prox_math.h"
#include <prox_update_inertia_tensor.h>

#include <tiny_is_number.h>
#include <tiny_is_finite.h>
#include <tiny_matrix_functions.h>

#include <cassert>

namespace prox
{

  // 2009-08-13 Kenny code reivew: Optimization replace diagonal6x6_type with diagonal_mass_type, maybe wait to optimize until all it working

template <typename T, typename Iterator>
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
}

template <typename T, typename Iterator>
inline void get_inverse_mass_matrix_eigen(Iterator begin, Iterator end,
                                          Eigen::SparseMatrix<T>& W)
{
    size_t const N = std::distance(begin, end);
    constexpr int block_size = 6;
    size_t total_size = N * block_size;

    std::vector<Eigen::Triplet<T>> triplets;
    triplets.reserve(N * block_size
                     * block_size); // Reserve for worst-case (all blocks full)

    size_t index = 0u;
    for (auto body = begin; body != end; ++body, ++index)
    {
        T inv_mass = 0;
        Eigen::Matrix3<T> inv_I = Eigen::Matrix3<T>::Zero();

        if (!body->is_fixed() && !body->is_scripted())
        {
            assert(std::abs(body->get_mass()) > 0
                   && "get_inverse_mass_matrix(): Divide by zero!");
            inv_mass = 1 / body->get_mass();

            // Check for valid inv_mass
            assert(std::isfinite(inv_mass)
                   && "get_inverse_mass_matrix(): Invalid mass");
            assert(inv_mass > 0 && "get_inverse_mass_matrix(): Negative mass");

            auto I_bf = body->get_inertia_bf();
            auto R = body->get_orientation().toRotationMatrix();

            // Compute inertia tensor in world frame and invert
            Eigen::Matrix3<T> I_wf = R * I_bf * R.transpose();
            inv_I = I_wf.inverse();
        }

        // Add the 6x6 block for the current body
        int base_row = index * block_size;

        // Mass terms (diagonal)
        for (int i = 0; i < 3; ++i)
        {
            triplets.emplace_back(base_row + i, base_row + i, inv_mass);
        }

        // Inertia tensor terms (3x3 sub-block)
        for (int i = 0; i < 3; ++i)
        {
            for (int j = 0; j < 3; ++j)
            {
                triplets.emplace_back(base_row + 3 + i, base_row + 3 + j,
                                      inv_I(i, j));
            }
        }
    }

    // Build the sparse matrix
    W.resize(total_size, total_size);
    W.setFromTriplets(triplets.begin(), triplets.end());
}
} // namespace prox
// PROX_GET_INVERSE_MASS_MATRIX_H
#endif
