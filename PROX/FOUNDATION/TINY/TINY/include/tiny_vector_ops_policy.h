#ifndef TINY_VECTOR_OPS_POLICY_H
#define TINY_VECTOR_OPS_POLICY_H

#include <cmath>

namespace tiny
{
  namespace detail
  {

    /**
     * VectorOpsPolicy class
     * This class implements the default policy specific for vectors.
     *
     * @tparam V    The vector type.
     */
    template <typename V>
    class VectorOpsPolicy
      {
      protected:
          using type_traits = typename V::type_traits;
          using value_traits = typename V::value_traits;

          enum
          {
              stride = type_traits::stride
          };
          enum
          {
              J = V::J
          };

      public:
          using real_type = typename V::real_type;
          using op_type = typename V::op_type;

          static real_type inner_prod(V const& lhs, V const& rhs)
          {
              op_type dot = type_traits::set_op_type(0);

              for (size_t j = 0; j < J; j += stride)
                  type_traits::add_assign(dot,
                                          type_traits::mul(V::accessor::cast(lhs, 0, j), V::accessor::cast(rhs, 0, j)));
              return type_traits::sum(dot);
        }

        static V cross (V const & lhs, V const & rhs)
        {
          real_type a,b,c;
          a = lhs[1]*rhs[2] - lhs[2]*rhs[1];
          b = lhs[2]*rhs[0] - lhs[0]*rhs[2];
          c = lhs[0]*rhs[1] - lhs[1]*rhs[0];
          return V::make(a,b,c);
        }

        static real_type norm (V const & m)
        {
          using std::sqrt;

          return sqrt( inner_prod(m,m) );
        }


        static V unit(V const & m)
        {
          V res(m);
          real_type const lgth     = norm(m);
          if(lgth>0)
          {
            op_type   const inv_norm = type_traits::set_op_type(1/lgth);
            for (size_t j = 0 ; j < J ; j+=stride )
              type_traits::mul_assign(V::accessor::cast(res,0,j),inv_norm);
          }
          return res;
        }

      };

  } //namespace detail
} //namespace tiny

// TINY_VECTOR_OPS_POLICY_H
#endif
