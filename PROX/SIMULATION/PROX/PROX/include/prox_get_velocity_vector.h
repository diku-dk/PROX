#ifndef PROX_GET_VELOCITY_VECTOR_H
#define PROX_GET_VELOCITY_VECTOR_H

#include "prox_math.h"
#include <tiny_is_number.h>

namespace prox
{

template <typename T, typename Iterator>
void get_velocity_vector(Iterator begin, Iterator end, NCVec6<T>& u)
{

    size_t const N = std::distance(begin,end);
    u.resize( N );

    size_t k = 0u;

    for (auto body = begin; body != end; ++body, ++k)
    {
        const auto& V = body->get_velocity();
        const auto& W = body->get_spin();

        assert(is_number(V(0)) || !"get_velocity_vector(): non number encountered");
        assert(is_number(V(1)) || !"get_velocity_vector(): non number encountered");
        assert(is_number(V(2)) || !"get_velocity_vector(): non number encountered");
        assert(is_number(W(0)) || !"get_velocity_vector(): non number encountered");
        assert(is_number(W(1)) || !"get_velocity_vector(): non number encountered");
        assert(is_number(W(2)) || !"get_velocity_vector(): non number encountered");

        auto& b = u(k);

        if (body->is_fixed())
        {
            b(0) = 0;
            b(1) = 0;
            b(2) = 0;
            b(3) = 0;
            b(4) = 0;
            b(5) = 0;
        }
      else
      {
        b(0) = V(0);
        b(1) = V(1);
        b(2) = V(2);
        b(3) = W(0);
        b(4) = W(1);
        b(5) = W(2);
      }
    }
}
} // namespace prox
// PROX_GET_VELOCITY_VECTOR_H
#endif
