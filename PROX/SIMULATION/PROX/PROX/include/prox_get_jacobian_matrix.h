#ifndef PROX_GET_JACOBIAN_MATRIX_H
#define PROX_GET_JACOBIAN_MATRIX_H

#include "prox_math.h"
#include <prox_property.h>
#include <prox_rigid_body.h>

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
        auto n_k = fromEigen(contact->normal);
        EigenVector3<T> t_k, s_k;

        size_t const mat_i = contact->bodyI->get_material_idx();
        size_t const mat_j = contact->bodyJ->get_material_idx();

        auto const* property = &(properties[mat_i][mat_j]);

        // 2012-07-15 Kenny code review: For box model it would make sense to have an option
        // that could pick s_t as the sliding direction vector! Not sure how to design prox
        // interfaces to accomodate this?

        if (property->is_isotropic())
        {
            tiny::orthonormal_vectors(
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

            s_k = tiny::rotate(Q, property->get_s_vector());
            float const c = tiny::inner_prod(n_k, s_k);
            s_k = tiny::unit(s_k - n_k * c);
            float const s_test = tiny::inner_prod(s_k, s_k);
            if (s_test > 10e-5f) { t_k = tiny::cross(n_k, s_k); }
            else
            {
                tiny::orthonormal_vectors(
                    s_k, t_k,
                    n_k); // 200X-YY-ZZ Kenny: TODO we might want to do this a little more clever?
            }
        }

        auto& J_ki = J(k, contact->bodyI->get_idx());

        // Fill in J_ki
        auto r_ki = contact->position - contact->bodyI->get_position();
        auto iXn = tiny::cross(r_ki, n_k);
        auto iXt = tiny::cross(r_ki, t_k);
        auto iXs = tiny::cross(r_ki, s_k);

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
        auto r_kj = contact->position - contact->bodyJ->get_position();
        auto jXn = tiny::cross(r_kj, n_k);
        auto jXt = tiny::cross(r_kj, t_k);
        auto jXs = tiny::cross(r_kj, s_k);

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

} // namespace prox

// PROX_GET_JACOBIAN_MATRIX_H
#endif
