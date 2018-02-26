#ifndef BIG_PROJECTED_GAUSS_SEIDEL_H
#define BIG_PROJECTED_GAUSS_SEIDEL_H

#include <big_types.h>
#include <big_return_codes.h>
#include <big_relative_convergence.h>
#include <big_absolute_convergence.h>
#include <big_prod_row.h>

#include <big_value_traits.h>
#include <big_is_number.h>

#include <cmath>
#include <stdexcept>
#include <cassert>

namespace big
{

      template < typename T>
      inline void projected_gauss_seidel(
        ublas::compressed_matrix<T> const & A
        , ublas::vector<T> const & b
        , boost::numeric::ublas::vector<T> & x
        , size_t const & max_iter
        , T      const & tol_abs
        , T      const & tol_rel
        , T      const & tol_stag
        , size_t       & status
        , size_t       & iter
        , T            & err
        , ublas::vector<T> * profiling = 0
        )  
      {
        using std::fabs;
        using std::min;
        using std::max;

        typedef          ublas::vector<T>                  vector_type;
        typedef typename vector_type::size_type            size_type;
        typedef          big::ValueTraits<T>               value_traits;

        if(max_iter <= 0u)
          throw std::invalid_argument("max_iterations must be larger than zero");

        if(tol_abs < value_traits::zero() )
          throw std::invalid_argument("absolute_tolerance must be non-negative");

        if(tol_rel < value_traits::zero() )
          throw std::invalid_argument("relative_tolerance must be non-negative");

        if(tol_stag < value_traits::zero() )
          throw std::invalid_argument("stagnation_tolerance must be non-negative");

        if(profiling == &x)
          throw std::logic_error("profiling must not point to x-vector");

        if(profiling == &b)
          throw std::logic_error("profiling must not point to b-vector");

        err      = value_traits::infinity();
        iter     = 0u;
        status   = OK;

        size_type const N = b.size();

        if(N==0)
          return;

        if(x.size()!=N)
          throw std::invalid_argument("size of x-vector were differnt from b?");

        if(profiling)
        {
          (*profiling).resize( max_iter );
          (*profiling).clear();
        }

        // Create a boolean flag indicating whether we need to compute
        // the value of the merit function or not.
        bool const compute_merit = profiling || (tol_abs > value_traits::zero()) || (tol_rel > value_traits::zero());

        status = ITERATING;

        for ( ;iter < max_iter;  ++iter)
        {
          T max_dx = value_traits::zero();

          for (size_type i = 0u; i < N; ++ i)
          {
            T y_i = b(i) + big::prod_row(A,x,i);

            T const A_ii = A(i,i);

            assert(A_ii> value_traits::zero() || A_ii<value_traits::zero() || !"projected_gauss_seidel: diagonal entry is zero?");

            T const old_x = x(i);
            T const x_tmp = - (y_i / A_ii) + old_x;

            assert(is_number(x_tmp) || !"projected_gauss_seidel: not a number encountered");

            T const x_i =  (x_tmp < value_traits::zero() ) ? value_traits::zero() :  x_tmp ;

            assert(is_number(x_i) || !"projected_gauss_seidel: not a number encountered");

            x(i) = x_i;

            max_dx = max( max_dx, fabs( x_i - old_x ) );
          }

          // Check for stagnation
          if(max_dx < tol_stag)
          {
            status = STAGNATION;
            return;
          }

          // Compute new merit value and perform check stopping criteria's
          if( compute_merit )
          {
            T const old_err = err;
            err = ublas::inner_prod(x,x)*value_traits::half();

            if(profiling)
              (*profiling)(iter) = err;

            if( absolute_convergence( err, tol_abs) )
            {
              status = ABSOLUTE_CONVERGENCE;
              return;
            }

            if( relative_convergence(old_err, err, tol_rel) )
            {
              status = RELATIVE_CONVERGENCE;
              return;
            }
          }
        }
      }

} // namespace big

// BIG_PROJECTED_GAUSS_SEIDEL_H
#endif
