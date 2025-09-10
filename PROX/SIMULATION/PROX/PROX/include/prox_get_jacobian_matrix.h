#ifndef PROX_GET_JACOBIAN_MATRIX_H
#define PROX_GET_JACOBIAN_MATRIX_H

#include "eigenhelperfunctions.h"
#include <prox_property.h>
#include <prox_rigid_body.h>

namespace prox
{
template <typename T, typename ContactIt, typename BodyContainer>
requires std::is_floating_point_v<T>
inline void get_jacobian_matrix_eigen(
    ContactIt begin, ContactIt end, const BodyContainer& bodies,
    const std::vector<std::vector<Property<T>>>& properties,
    Eigen::SparseMatrix<T>& J, size_t K)
{
    size_t const N = bodies.size();

    // Resize the sparse matrix
    J.resize(4 * K, 6 * N); // 4 constraints per contact, 6 DOFs per body

    std::vector<Eigen::Triplet<T>>
        triplets; // To store non-zero entries in the sparse matrix

    size_t k = 0u;
    for (auto contact = begin; contact != end; ++contact, ++k)
    {
        auto n_k = contact->normal;
        EigenVector3<T> t_k, s_k;

        size_t const mat_i = contact->bodyI->get_material_idx();
        size_t const mat_j = contact->bodyJ->get_material_idx();

        auto const* property = &(properties[mat_i][mat_j]);

        if (property->is_isotropic()) { orthonormal_vectors(s_k, t_k, n_k); }
        else
        {
            assert(false || !"Anisotropic friction not tested");
            EigenQuaternion<T> Q;
            size_t const master = property->get_master_material_idx();
            if (contact->bodyI->get_material_idx() == master)
            {
                Q = contact->bodyI->get_orientation();
            }
            else if (contact->bodyJ->get_material_idx() == master)
            {
                Q = contact->bodyJ->get_orientation();
            }
            else { assert(false || !"Bad master idx"); }

            s_k = rotate(Q, property->get_s_vector());
            float const c = dot(n_k, s_k);
            s_k = (s_k - n_k * c).normalized();
            float const s_test = dot(s_k, s_k);
            if (s_test > 10e-5f) { t_k = cross(n_k, s_k); }
            else { orthonormal_vectors(s_k, t_k, n_k); }
        }

        // Access the entries for body I
        size_t i_idx = contact->bodyI->get_idx();
        auto& J_ki = J; // Sparse matrix reference
        EigenVector3<T> r_ki
            = contact->position - contact->bodyI->get_position();
        auto iXn = cross(r_ki, n_k);
        auto iXt = cross(r_ki, t_k);
        auto iXs = cross(r_ki, s_k);

        // Fill in J_ki
        triplets.emplace_back(4 * k + 0, 6 * i_idx + 0, -n_k(0));
        triplets.emplace_back(4 * k + 0, 6 * i_idx + 1, -n_k(1));
        triplets.emplace_back(4 * k + 0, 6 * i_idx + 2, -n_k(2));
        triplets.emplace_back(4 * k + 0, 6 * i_idx + 3, -iXn(0));
        triplets.emplace_back(4 * k + 0, 6 * i_idx + 4, -iXn(1));
        triplets.emplace_back(4 * k + 0, 6 * i_idx + 5, -iXn(2));

        triplets.emplace_back(4 * k + 1, 6 * i_idx + 0, -t_k(0));
        triplets.emplace_back(4 * k + 1, 6 * i_idx + 1, -t_k(1));
        triplets.emplace_back(4 * k + 1, 6 * i_idx + 2, -t_k(2));
        triplets.emplace_back(4 * k + 1, 6 * i_idx + 3, -iXt(0));
        triplets.emplace_back(4 * k + 1, 6 * i_idx + 4, -iXt(1));
        triplets.emplace_back(4 * k + 1, 6 * i_idx + 5, -iXt(2));

        triplets.emplace_back(4 * k + 2, 6 * i_idx + 0, -s_k(0));
        triplets.emplace_back(4 * k + 2, 6 * i_idx + 1, -s_k(1));
        triplets.emplace_back(4 * k + 2, 6 * i_idx + 2, -s_k(2));
        triplets.emplace_back(4 * k + 2, 6 * i_idx + 3, -iXs(0));
        triplets.emplace_back(4 * k + 2, 6 * i_idx + 4, -iXs(1));
        triplets.emplace_back(4 * k + 2, 6 * i_idx + 5, -iXs(2));

        triplets.emplace_back(4 * k + 3, 6 * i_idx + 0, (0.0f));
        triplets.emplace_back(4 * k + 3, 6 * i_idx + 1, (0.0f));
        triplets.emplace_back(4 * k + 3, 6 * i_idx + 2, (0.0f));
        triplets.emplace_back(4 * k + 3, 6 * i_idx + 3, -n_k(0));
        triplets.emplace_back(4 * k + 3, 6 * i_idx + 4, -n_k(1));
        triplets.emplace_back(4 * k + 3, 6 * i_idx + 5, -n_k(2));

        // Access the entries for body J (same as for body I)
        size_t j_idx = contact->bodyJ->get_idx();
        auto& J_kj = J;
        EigenVector3<T> r_kj
            = contact->position - contact->bodyJ->get_position();
        auto jXn = cross(r_kj, n_k);
        auto jXt = cross(r_kj, t_k);
        auto jXs = cross(r_kj, s_k);

        // Fill in J_kj
        triplets.emplace_back(4 * k + 0, 6 * j_idx + 0, n_k(0));
        triplets.emplace_back(4 * k + 0, 6 * j_idx + 1, n_k(1));
        triplets.emplace_back(4 * k + 0, 6 * j_idx + 2, n_k(2));
        triplets.emplace_back(4 * k + 0, 6 * j_idx + 3, jXn(0));
        triplets.emplace_back(4 * k + 0, 6 * j_idx + 4, jXn(1));
        triplets.emplace_back(4 * k + 0, 6 * j_idx + 5, jXn(2));

        triplets.emplace_back(4 * k + 1, 6 * j_idx + 0, t_k(0));
        triplets.emplace_back(4 * k + 1, 6 * j_idx + 1, t_k(1));
        triplets.emplace_back(4 * k + 1, 6 * j_idx + 2, t_k(2));
        triplets.emplace_back(4 * k + 1, 6 * j_idx + 3, jXt(0));
        triplets.emplace_back(4 * k + 1, 6 * j_idx + 4, jXt(1));
        triplets.emplace_back(4 * k + 1, 6 * j_idx + 5, jXt(2));

        triplets.emplace_back(4 * k + 2, 6 * j_idx + 0, s_k(0));
        triplets.emplace_back(4 * k + 2, 6 * j_idx + 1, s_k(1));
        triplets.emplace_back(4 * k + 2, 6 * j_idx + 2, s_k(2));
        triplets.emplace_back(4 * k + 2, 6 * j_idx + 3, jXs(0));
        triplets.emplace_back(4 * k + 2, 6 * j_idx + 4, jXs(1));
        triplets.emplace_back(4 * k + 2, 6 * j_idx + 5, jXs(2));

        triplets.emplace_back(4 * k + 3, 6 * j_idx + 0, T(0.0f));
        triplets.emplace_back(4 * k + 3, 6 * j_idx + 1, T(0.0f));
        triplets.emplace_back(4 * k + 3, 6 * j_idx + 2, T(0.0f));
        triplets.emplace_back(4 * k + 3, 6 * j_idx + 3, n_k(0));
        triplets.emplace_back(4 * k + 3, 6 * j_idx + 4, n_k(1));
        triplets.emplace_back(4 * k + 3, 6 * j_idx + 5, n_k(2));
    }

    // Once all triplets are populated, set them in the sparse matrix
    J.setFromTriplets(triplets.begin(), triplets.end());
}

} // namespace prox

// PROX_GET_JACOBIAN_MATRIX_H
#endif
