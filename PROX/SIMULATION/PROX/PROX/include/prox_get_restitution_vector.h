#ifndef RPOX_GET_RESTITUTION_VECTOR_H
#define RPOX_GET_RESTITUTION_VECTOR_H

#include "prox_math.h"
#include "prox_property.h"
#include <util_log.h>

#include <cstddef>

namespace prox
{

template <typename T, typename Iterator>
void get_restitution_vector(
    Iterator begin, Iterator end,
    const std::vector<std::vector<Property<T>>>& properties, NCVec4<T>& e,
    size_t K)
{
    util::Log logging;

    e.resize(K);

    logging << "get_restitution_vector(): bounce is on" << util::Log::newline();

    size_t k = 0u;
    for (auto contact = begin; contact != end; ++contact, ++k)
    {
        auto material_i = contact->bodyI->get_material_idx();
        auto material_j = contact->bodyJ->get_material_idx();

        auto e_k
            = properties[material_i][material_j].get_restitution_coefficient();

        auto& b = e(k);

        b(0) = e_k;
        b(1) = 0;
        b(2) = 0;
        b(3) = 0;
    }
}

} // namespace prox

// RPOX_GET_RESTITUTION_VECTOR_H
#endif
