#ifndef PROX_GET_INVERSE_MASS_MATRIX_H
#define PROX_GET_INVERSE_MASS_MATRIX_H

#include "eigenhelperfunctions.h"
#include <prox_update_inertia_tensor.h>

#include <eigenhelperfunctions.h>

#include <cassert>

namespace prox
{

  // 2009-08-13 Kenny code reivew: Optimization replace diagonal6x6_type with diagonal_mass_type, maybe wait to optimize until all it working

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

            inv_mass = 1 / (body->get_mass() /** 1.5*/);

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

} // namespace prox
// PROX_GET_INVERSE_MASS_MATRIX_H
#endif
