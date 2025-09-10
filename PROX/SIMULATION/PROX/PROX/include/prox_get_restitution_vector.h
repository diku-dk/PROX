#ifndef RPOX_GET_RESTITUTION_VECTOR_H
#define RPOX_GET_RESTITUTION_VECTOR_H

#include "eigenhelperfunctions.h"
#include "prox_property.h"
#include <util_log.h>

#include <cstddef>

namespace prox
{


template <typename T, typename Iterator>
void get_restitution_vector_eigen(
    Iterator begin, Iterator end,
    const std::vector<std::vector<Property<T>>>& properties,
    Eigen::VectorX<T>& e, size_t K)
{
    util::Log logging;

    e.resize(K * 4);  // Resize to (K * 4) x 1
    e.setZero();      // Initialize all elements to 0

    logging << "get_restitution_vector(): bounce is on" << util::Log::newline();

    size_t k = 0;
    for (auto contact = begin; contact != end; ++contact, ++k)
    {
        auto material_i = contact->bodyI->get_material_idx();
        auto material_j = contact->bodyJ->get_material_idx();

        auto e_k
            = properties[material_i][material_j].get_restitution_coefficient();

        // Calculate the starting index for this contact's data
        size_t idx = k * 4;

        // Set the first element to e_k, others remain 0
        e(idx) = e_k;
    }
}

} // namespace prox

// RPOX_GET_RESTITUTION_VECTOR_H
#endif
