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

template <typename T, typename Iterator>
void get_velocity_vector_eigen(Iterator begin, Iterator end,
                               Eigen::VectorX<T>& q)
{
    size_t N = std::distance(begin, end);
    q.resize(N * 6);  // Resize to (N*6)x1

    size_t k = 0;
    for (auto body = begin; body != end; ++body, ++k)
    {
        const auto& V = body->get_velocity();
        const auto& W = body->get_spin();

        assert(is_number(V(0))
               || !"get_velocity_vector(): non number encountered");
        assert(is_number(V(1))
               || !"get_velocity_vector(): non number encountered");
        assert(is_number(V(2))
               || !"get_velocity_vector(): non number encountered");
        assert(is_number(W(0))
               || !"get_velocity_vector(): non number encountered");
        assert(is_number(W(1))
               || !"get_velocity_vector(): non number encountered");
        assert(is_number(W(2))
               || !"get_velocity_vector(): non number encountered");

        // Calculate the starting index for this body's data
        size_t idx = k * 6;
        if (body->is_fixed())
        {
            q(idx) = 0;      // x-position
            q(idx + 1) = 0;      // y-position
            q(idx + 2) = 0;      // z-position
            q(idx + 3) = 0;     // Quaternion w
            q(idx + 4) = 0;     // Quaternion x
            q(idx + 5) = 0;     // Quaternion y
        }
        else
        {
            q(idx) = V(0);      // x-position
            q(idx + 1) = V(1);      // y-position
            q(idx + 2) = V(2);      // z-position
            q(idx + 3) = W(0);     // Quaternion w
            q(idx + 4) = W(1);     // Quaternion x
            q(idx + 5) = W(2);     // Quaternion y
        }
    }
}
} // namespace prox
// PROX_GET_VELOCITY_VECTOR_H
#endif
