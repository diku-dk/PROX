#ifndef PROX_GET_INVERSE_MASS_MATRIX_H
#define PROX_GET_INVERSE_MASS_MATRIX_H

#include "prox_math.h"
#include <prox_update_inertia_tensor.h>

#include <tiny_is_number.h>
#include <tiny_is_finite.h>
#include <tiny_matrix_functions.h>

#include <cassert>

namespace prox
{

  // 2009-08-13 Kenny code reivew: Optimization replace diagonal6x6_type with diagonal_mass_type, maybe wait to optimize until all it working

template <typename T, typename Iterator>
inline void get_inverse_mass_matrix(Iterator begin, Iterator end, DiagonalMatrix<6, T>& W)
{
    size_t const N = std::distance(begin,end);

    W.resize( N );

    size_t index = 0u;
    for (auto body = begin; body != end; ++body, ++index)
    {
        T inv_mass = 0;
        EigenMatrix3<T> inv_I{
            {0, 0, 0},
            {0, 0, 0},
            {0, 0, 0}
        };

        if (!body->is_fixed() && !body->is_scripted())
        {
            assert(fabs(body->get_mass()) > 0 || !"get_inverse_mass_matrix(): Divide by zero!");

            inv_mass = 1 / body->get_mass();

            assert(is_number(inv_mass) || !"get_inverse_mass_matrix(): Nan");
            assert(is_finite(inv_mass) || !"get_inverse_mass_matrix(): Inf");
            assert(inv_mass > 0 || !"get_inverse_mass_matrix(): Negative mass");

            auto I_bf = body->get_inertia_bf();
            auto R = EigenMatrix3<T>(body->get_orientation());

            detail::update_inertia_tensor(R, I_bf, inv_I);

            inv_I = inv_I.inverse().eval();
      }

      auto& b = W(index);

      b(0,0) = inv_mass;
      b(1,1) = inv_mass;
      b(2,2) = inv_mass;
      b(3,3) = inv_I(0,0);
      b(3,4) = inv_I(0,1);
      b(3,5) = inv_I(0,2);
      b(4,3) = inv_I(1,0);
      b(4,4) = inv_I(1,1);
      b(4,5) = inv_I(1,2);
      b(5,3) = inv_I(2,0);
      b(5,4) = inv_I(2,1);
      b(5,5) = inv_I(2,2);
    }
}
} // namespace prox
// PROX_GET_INVERSE_MASS_MATRIX_H
#endif
