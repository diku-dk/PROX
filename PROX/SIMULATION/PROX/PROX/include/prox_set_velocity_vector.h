#ifndef PROX_SET_VELOCITY_VECTOR_H
#define PROX_SET_VELOCITY_VECTOR_H

#include "prox_math.h"
#include "tiny_math_types.h"
#include <tiny_is_number.h>

namespace prox
{

template <typename T, typename Iterator>
inline void set_velocity_vector(Iterator begin, Iterator end,
                                const NCVec6<T>& u)
{
    EigenVector3<T> V, W;
    size_t const N = std::distance(begin, end);

    assert(u.size() == N
           || !"set_velocity_vector(): u has incorrect dimension");

    size_t k = 0u;
    for (auto body = begin; body != end; ++body, ++k)
    {
        if (body->is_fixed() || body->is_scripted()) continue;

        auto const& b = u(k);

        V(0) = b(0);
        V(1) = b(1);
        V(2) = b(2);

        W(0) = b(3);
        W(1) = b(4);
        W(2) = b(5);

        assert(is_number(V(0))
               || !"set_velocity_vector(): non number encountered");
        assert(is_number(V(1))
               || !"set_velocity_vector(): non number encountered");
        assert(is_number(V(2))
               || !"set_velocity_vector(): non number encountered");
        assert(is_number(W(0))
               || !"set_velocity_vector(): non number encountered");
        assert(is_number(W(1))
               || !"set_velocity_vector(): non number encountered");
        assert(is_number(W(2))
               || !"set_velocity_vector(): non number encountered");

        body->set_velocity(V);
        body->set_spin(W);
    }
}
} // namespace prox

// PROX_SET_VELOCITY_VECTOR_H
#endif
