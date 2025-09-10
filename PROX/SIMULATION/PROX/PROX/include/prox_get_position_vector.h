#ifndef PROX_GET_POSITION_VECTOR_H
#define PROX_GET_POSITION_VECTOR_H

#include "eigenhelperfunctions.h"

namespace prox
{

template <typename T, typename Iterator>
void get_position_vector_eigen(Iterator begin, Iterator end,
                               Eigen::VectorX<T>& q)
{
    size_t N = std::distance(begin, end);
    q.resize(N * 7);  // Resize to (N*7)x1

    size_t k = 0;
    for (auto body = begin; body != end; ++body, ++k)
    {
        const auto& r = body->get_position();
        const auto& Q = body->get_orientation();

        assert(is_number(r(0))
               || !"get_position_vector(): non number encountered");
        assert(is_number(r(1))
               || !"get_position_vector(): non number encountered");
        assert(is_number(r(2))
               || !"get_position_vector(): non number encountered");

        // Calculate the starting index for this body's data
        size_t idx = k * 7;
        q(idx) = r(0);      // x-position
        q(idx + 1) = r(1);      // y-position
        q(idx + 2) = r(2);      // z-position
        q(idx + 3) = Q.w();     // Quaternion w
        q(idx + 4) = Q.x();     // Quaternion x
        q(idx + 5) = Q.y();     // Quaternion y
        q(idx + 6) = Q.z();     // Quaternion z
    }
}
}// namespace prox
// PROX_GET_POSITION_VECTOR_H
#endif
