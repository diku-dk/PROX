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
} // namespace prox
// PROX_SET_POSITION_VECTOR_H
#endif
