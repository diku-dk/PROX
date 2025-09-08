#ifndef PROX_SET_POSITION_VECTOR_H
#define PROX_SET_POSITION_VECTOR_H

#include "prox_math.h"
#include "tiny_math_types.h"
#include <tiny_is_number.h>

namespace prox
{

template <typename T, typename Iterator>
inline void set_position_vector(Iterator begin, Iterator end, NCVec7<T>& q)
{
    size_t const N = std::distance(begin, end);

    assert(q.size() == N || !"set_position_vector(): q has incorrect dimension");

    size_t k = 0u;
    for (auto body = begin; body != end; ++body, ++k)
    {
      if( body->is_fixed() || body->is_scripted() )
        continue;

      EigenVector3<T> r;
      EigenQuaternion<T> Q;

      const auto& b = q(k);

      r(0) = b(0);
      r(1) = b(1);
      r(2) = b(2);
      Q.w() = b(3);
      Q.x() = b(4);
      Q.y() = b(5);
      Q.z() = b(6);

      assert(is_number(r(0))        || !"set_position_vector(): non number encountered");
      assert(is_number(r(1))        || !"set_position_vector(): non number encountered");
      assert(is_number(r(2))        || !"set_position_vector(): non number encountered");
      body->set_position( r );
      body->set_orientation( Q );
    }
}

template <typename T, typename Iterator>
inline void set_position_vector_eigen(Iterator begin, Iterator end,
                                      Eigen::VectorX<T>& q)
{
    size_t const N = std::distance(begin, end);

    // Check if q has the correct dimension (N*7)
    assert(q.size() == N * 7
           && "set_position_vector(): q has incorrect dimension");

    size_t k = 0;
    for (auto body = begin; body != end; ++body, ++k)
    {
        if (body->is_fixed() || body->is_scripted()) continue;

        // Calculate the starting index for this body's data
        size_t idx = k * 7;

        // Extract position and quaternion from the flat vector
        Eigen::Matrix<T, 3, 1> r;
        r(0) = q(idx);
        r(1) = q(idx + 1);
        r(2) = q(idx + 2);

        Eigen::Quaternion<T> Q;
        Q.w() = q(idx + 3);
        Q.x() = q(idx + 4);
        Q.y() = q(idx + 5);
        Q.z() = q(idx + 6);

        // Check for valid numbers
        assert(is_number(r(0))
               && "set_position_vector(): non number encountered");
        assert(is_number(r(1))
               && "set_position_vector(): non number encountered");
        assert(is_number(r(2))
               && "set_position_vector(): non number encountered");
        assert(is_number(Q.w())
               && "set_position_vector(): non number encountered");
        assert(is_number(Q.x())
               && "set_position_vector(): non number encountered");
        assert(is_number(Q.y())
               && "set_position_vector(): non number encountered");
        assert(is_number(Q.z())
               && "set_position_vector(): non number encountered");

        body->set_position(r);
        body->set_orientation(Q);
    }
}
} // namespace prox
// PROX_SET_POSITION_VECTOR_H
#endif
