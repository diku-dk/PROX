#ifndef PROX_GET_JACOBIAN_MATRIX_H
#define PROX_GET_JACOBIAN_MATRIX_H

#include "prox_math.h"
#include <prox_property.h>
#include <prox_rigid_body.h>
#include <eigen3/Eigen/Sparse>

namespace prox
{

template <typename T, typename ContactIt, typename BodyContainer>
requires std::is_floating_point_v<T>
inline void get_jacobian_matrix(ContactIt begin, ContactIt end, const BodyContainer& bodies,
                                const std::vector<std::vector<Property<T>>>& properties,
                                CRMatrix<4, 6, T>& J, size_t K)
{
    size_t const N = bodies.size();

    J.resize(K, N, 2 * K);

    size_t k = 0u;

    for (auto contact = begin; contact != end; ++contact, ++k)
    {
        auto n_k = contact->normal;
        EigenVector3<T> t_k, s_k;

        size_t const mat_i = contact->bodyI->get_material_idx();
        size_t const mat_j = contact->bodyJ->get_material_idx();

        auto const* property = &(properties[mat_i][mat_j]);

        // 2012-07-15 Kenny code review: For box model it would make sense to have an option
        // that could pick s_t as the sliding direction vector! Not sure how to design prox
        // interfaces to accomodate this?

        if (property->is_isotropic())
        {
            orthonormal_vectors(
                s_k, t_k,
                n_k); // 200X-YY-ZZ Kenny: TODO we might want to do this a little more clever?
        }
        else
        {
            assert(false || !"get_jacobian_matrix(): Anisotropic friction is not yet tested");

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
            else { assert(false || !"get_jacobian_matrix(): master idx was bad"); }

            s_k = rotate(Q, property->get_s_vector());
            float const c = dot(n_k, s_k);
            s_k = (s_k - n_k * c).normalized();
            float const s_test = dot(s_k, s_k);
            if (s_test > 10e-5f) { t_k = cross(n_k, s_k); }
            else
            {
                orthonormal_vectors(
                    s_k, t_k,
                    n_k); // 200X-YY-ZZ Kenny: TODO we might want to do this a little more clever?
            }
        }

        auto& J_ki = J(k, contact->bodyI->get_idx());

        // Fill in J_ki
        EigenVector3<T> r_ki
            = contact->position - contact->bodyI->get_position();
        auto iXn = cross(r_ki, n_k);
        auto iXt = cross(r_ki, t_k);
        auto iXs = cross(r_ki, s_k);

        J_ki(0, 0) = -n_k(0);
        J_ki(0, 1) = -n_k(1);
        J_ki(0, 2) = -n_k(2);
        J_ki(0, 3) = -iXn(0);
        J_ki(0, 4) = -iXn(1);
        J_ki(0, 5) = -iXn(2);
        J_ki(1, 0) = -t_k(0);
        J_ki(1, 1) = -t_k(1);
        J_ki(1, 2) = -t_k(2);
        J_ki(1, 3) = -iXt(0);
        J_ki(1, 4) = -iXt(1);
        J_ki(1, 5) = -iXt(2);
        J_ki(2, 0) = -s_k(0);
        J_ki(2, 1) = -s_k(1);
        J_ki(2, 2) = -s_k(2);
        J_ki(2, 3) = -iXs(0);
        J_ki(2, 4) = -iXs(1);
        J_ki(2, 5) = -iXs(2);
        J_ki(3, 0) = 0.0f;
        J_ki(3, 1) = 0.0f;
        J_ki(3, 2) = 0.0f;
        J_ki(3, 3) = -n_k(0);
        J_ki(3, 4) = -n_k(1);
        J_ki(3, 5) = -n_k(2);

        auto& J_kj = J(k, contact->bodyJ->get_idx());

        // Fill in J_kj
        EigenVector3<T> r_kj
            = contact->position - contact->bodyJ->get_position();
        auto jXn = cross(r_kj, n_k);
        auto jXt = cross(r_kj, t_k);
        auto jXs = cross(r_kj, s_k);

        J_kj(0, 0) = n_k(0);
        J_kj(0, 1) = n_k(1);
        J_kj(0, 2) = n_k(2);
        J_kj(0, 3) = jXn(0);
        J_kj(0, 4) = jXn(1);
        J_kj(0, 5) = jXn(2);
        J_kj(1, 0) = t_k(0);
        J_kj(1, 1) = t_k(1);
        J_kj(1, 2) = t_k(2);
        J_kj(1, 3) = jXt(0);
        J_kj(1, 4) = jXt(1);
        J_kj(1, 5) = jXt(2);
        J_kj(2, 0) = s_k(0);
        J_kj(2, 1) = s_k(1);
        J_kj(2, 2) = s_k(2);
        J_kj(2, 3) = jXs(0);
        J_kj(2, 4) = jXs(1);
        J_kj(2, 5) = jXs(2);
        J_kj(3, 0) = 0.0f;
        J_kj(3, 1) = 0.0f;
        J_kj(3, 2) = 0.0f;
        J_kj(3, 3) = n_k(0);
        J_kj(3, 4) = n_k(1);
        J_kj(3, 5) = n_k(2);
    }
}

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

        triplets.emplace_back(4 * k + 3, 6 * i_idx + 0, 0.0f);
        triplets.emplace_back(4 * k + 3, 6 * i_idx + 1, 0.0f);
        triplets.emplace_back(4 * k + 3, 6 * i_idx + 2, 0.0f);
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

        triplets.emplace_back(4 * k + 3, 6 * j_idx + 0, 0.0f);
        triplets.emplace_back(4 * k + 3, 6 * j_idx + 1, 0.0f);
        triplets.emplace_back(4 * k + 3, 6 * j_idx + 2, 0.0f);
        triplets.emplace_back(4 * k + 3, 6 * j_idx + 3, n_k(0));
        triplets.emplace_back(4 * k + 3, 6 * j_idx + 4, n_k(1));
        triplets.emplace_back(4 * k + 3, 6 * j_idx + 5, n_k(2));
    }

    // Once all triplets are populated, set them in the sparse matrix
    J.setFromTriplets(triplets.begin(), triplets.end());
}

template <typename T> using Mat4x6 = Eigen::Matrix<T, 4, 6>;
template <typename T> using Vec3 = Eigen::Matrix<T, 3, 1>;

template <typename T, typename ContactIt, typename BodyContainer>
requires std::is_floating_point_v<T>
inline void get_jacobian_matrix_sparse(
    ContactIt begin, ContactIt end, const BodyContainer& bodies,
    const std::vector<std::vector<Property<T>>>& properties,
    Eigen::SparseMatrix<T>& J, size_t K)
{
    static_assert(std::is_floating_point_v<T>, "T must be floating point");

    const size_t N = bodies.size();
    const Eigen::Index rows = static_cast<Eigen::Index>(4 * K);
    const Eigen::Index cols = static_cast<Eigen::Index>(6 * N);

    // Reserve triplets: each contact contributes up to 2 * (4*6) = 48 entries.
    std::vector<Eigen::Triplet<T>> triplets;
    triplets.reserve(std::distance(begin, end) * 48);

    size_t k = 0u;
    for (auto contact = begin; contact != end; ++contact, ++k)
    {
        Vec3<T> n_k = contact->normal;
        Vec3<T> t_k, s_k;

        size_t const mat_i = contact->bodyI->get_material_idx();
        size_t const mat_j = contact->bodyJ->get_material_idx();
        auto const* property = &(properties[mat_i][mat_j]);

        if (property->is_isotropic())
        {
            // provide or implement this helper; it must fill s_k,t_k orthonormal to n_k
            orthonormal_vectors(s_k, t_k, n_k);
        }
        else
        {
            assert(false
                   || !"get_jacobian_matrix(): Anisotropic friction is not yet "
                       "tested");
            Eigen::Quaternion<T> Q;
            size_t const master = property->get_master_material_idx();
            if (contact->bodyI->get_material_idx() == master)
            {
                Q = contact->bodyI->get_orientation();
            }
            else if (contact->bodyJ->get_material_idx() == master)
            {
                Q = contact->bodyJ->get_orientation();
            }
            else
            {
                assert(false || !"get_jacobian_matrix(): master idx was bad");
            }

            s_k = Q * property->get_s_vector();
            T c = n_k.dot(s_k);
            s_k = (s_k - n_k * c).normalized();
            T s_test = s_k.dot(s_k);
            if (s_test > T(1e-4)) { t_k = n_k.cross(s_k); }
            else { orthonormal_vectors(s_k, t_k, n_k); }
        }

        // Build the 4x6 block for body I
        Mat4x6<T> blockI;
        blockI.setZero();
        Vec3<T> r_ki = contact->position - contact->bodyI->get_position();
        Vec3<T> iXn = r_ki.cross(n_k);
        Vec3<T> iXt = r_ki.cross(t_k);
        Vec3<T> iXs = r_ki.cross(s_k);

        blockI(0, 0) = -n_k(0);
        blockI(0, 1) = -n_k(1);
        blockI(0, 2) = -n_k(2);
        blockI(0, 3) = -iXn(0);
        blockI(0, 4) = -iXn(1);
        blockI(0, 5) = -iXn(2);

        blockI(1, 0) = -t_k(0);
        blockI(1, 1) = -t_k(1);
        blockI(1, 2) = -t_k(2);
        blockI(1, 3) = -iXt(0);
        blockI(1, 4) = -iXt(1);
        blockI(1, 5) = -iXt(2);

        blockI(2, 0) = -s_k(0);
        blockI(2, 1) = -s_k(1);
        blockI(2, 2) = -s_k(2);
        blockI(2, 3) = -iXs(0);
        blockI(2, 4) = -iXs(1);
        blockI(2, 5) = -iXs(2);

        blockI(3, 0) = T(0);
        blockI(3, 1) = T(0);
        blockI(3, 2) = T(0);
        blockI(3, 3) = -n_k(0);
        blockI(3, 4) = -n_k(1);
        blockI(3, 5) = -n_k(2);

        // Build the 4x6 block for body J
        Mat4x6<T> blockJ;
        blockJ.setZero();
        Vec3<T> r_kj = contact->position - contact->bodyJ->get_position();
        Vec3<T> jXn = r_kj.cross(n_k);
        Vec3<T> jXt = r_kj.cross(t_k);
        Vec3<T> jXs = r_kj.cross(s_k);

        blockJ(0, 0) = n_k(0);
        blockJ(0, 1) = n_k(1);
        blockJ(0, 2) = n_k(2);
        blockJ(0, 3) = jXn(0);
        blockJ(0, 4) = jXn(1);
        blockJ(0, 5) = jXn(2);

        blockJ(1, 0) = t_k(0);
        blockJ(1, 1) = t_k(1);
        blockJ(1, 2) = t_k(2);
        blockJ(1, 3) = jXt(0);
        blockJ(1, 4) = jXt(1);
        blockJ(1, 5) = jXt(2);

        blockJ(2, 0) = s_k(0);
        blockJ(2, 1) = s_k(1);
        blockJ(2, 2) = s_k(2);
        blockJ(2, 3) = jXs(0);
        blockJ(2, 4) = jXs(1);
        blockJ(2, 5) = jXs(2);

        blockJ(3, 0) = T(0);
        blockJ(3, 1) = T(0);
        blockJ(3, 2) = T(0);
        blockJ(3, 3) = n_k(0);
        blockJ(3, 4) = n_k(1);
        blockJ(3, 5) = n_k(2);

        // Insert blockI into triplets: rows = 4*k..4*k+3, cols = 6*idxI..6*idxI+5
        const size_t idxI = contact->bodyI->get_idx();
        const size_t idxJ = contact->bodyJ->get_idx();
        const int row0 = static_cast<int>(4 * k);
        const int colI0 = static_cast<int>(6 * idxI);
        const int colJ0 = static_cast<int>(6 * idxJ);

        for (int r = 0; r < 4; ++r)
            for (int c = 0; c < 6; ++c)
                triplets.emplace_back(row0 + r, colI0 + c, blockI(r, c));

        for (int r = 0; r < 4; ++r)
            for (int c = 0; c < 6; ++c)
                triplets.emplace_back(row0 + r, colJ0 + c, blockJ(r, c));
    }

    J.resize(rows, cols);
    J.setFromTriplets(triplets.begin(), triplets.end());
    // Optionally compress
    J.makeCompressed();
}

} // namespace prox

// PROX_GET_JACOBIAN_MATRIX_H
#endif
