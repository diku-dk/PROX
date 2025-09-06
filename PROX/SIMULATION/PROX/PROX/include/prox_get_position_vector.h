#ifndef PROX_GET_POSITION_VECTOR_H
#define PROX_GET_POSITION_VECTOR_H

#include "prox_math.h"
#include <tiny_is_number.h>

namespace prox
{

template <typename T, typename Iterator>
void get_position_vector(Iterator begin, Iterator end, NCVec7<T>& q)
{
    size_t N = std::distance(begin,end);
    q.resize( N );

    size_t k = 0u;
    for (auto body = begin; body != end; ++body, ++k)
    {
        const auto& r = body->get_position();
        const auto& Q = body->get_orientation();

        assert(is_number(r(0)) || !"get_position_vector(): non number encountered");
        assert(is_number(r(1)) || !"get_position_vector(): non number encountered");
        assert(is_number(r(2)) || !"get_position_vector(): non number encountered");

        auto& b = q(k);

        b(0) = r(0);
        b(1) = r(1);
        b(2) = r(2);
        b(3) = Q.w();
        b(4) = Q.x();
        b(5) = Q.y();
        b(6) = Q.z();
    }
}
}// namespace prox
// PROX_GET_POSITION_VECTOR_H
#endif
