#ifndef PROX_GET_FRICTION_COEFFICIENT_VECTOR_H
#define PROX_GET_FRICTION_COEFFICIENT_VECTOR_H

#include "prox_math.h"
#include "prox_property.h"
#include <cstddef>

namespace prox
{

template <typename T, typename Iterator>
void get_friction_coefficient_vector(
    Iterator begin, Iterator end,
    const std::vector<std::vector<Property<T>>>& properties, NCVec4<T>& mu,
    size_t K)
{
    mu.resize(K);

    size_t k = 0u;
    for (auto contact = begin; contact != end; ++contact, ++k)
    {
        auto material_i = contact->bodyI->get_material_idx();
        auto material_j = contact->bodyJ->get_material_idx();

        auto const mu_s
            = properties[material_i][material_j].get_friction_coefficients()(0);
        auto const mu_t
            = properties[material_i][material_j].get_friction_coefficients()(1);
        auto const mu_tau
            = properties[material_i][material_j].get_friction_coefficients()(2);

        auto& b = mu(k);

        b(0) = 0;
        b(1) = mu_s;
        b(2) = mu_t;
        b(3) = mu_tau;
    }
}

template <typename T, typename Iterator>
void get_friction_coefficient_vector_eigen(
    Iterator begin, Iterator end,
    const std::vector<std::vector<Property<T>>>& properties,
    Eigen::VectorX<T>& mu, size_t K)
{
    mu.resize(K * 4);  // Resize to (K * 4) x 1

    size_t k = 0;
    for (auto contact = begin; contact != end; ++contact, ++k)
    {
        auto material_i = contact->bodyI->get_material_idx();
        auto material_j = contact->bodyJ->get_material_idx();

        auto const mu_s
            = properties[material_i][material_j].get_friction_coefficients()(0);
        auto const mu_t
            = properties[material_i][material_j].get_friction_coefficients()(1);
        auto const mu_tau
            = properties[material_i][material_j].get_friction_coefficients()(2);

        // Calculate the starting index for this contact's data
        size_t idx = k * 4;

        // Set the friction coefficients
        mu(idx) = 0;      // First element is 0
        mu(idx + 1) = mu_s;   // Second element is mu_s
        mu(idx + 2) = mu_t;   // Third element is mu_t
        mu(idx + 3) = mu_tau; // Fourth element is mu_tau
    }
}

} // namespace prox

// PROX_GET_FRICTION_COEFFICIENT_VECTOR_H
#endif
