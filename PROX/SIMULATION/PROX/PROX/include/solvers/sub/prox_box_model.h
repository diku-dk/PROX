#ifndef PROX_BOX_MODEL_H
#define PROX_BOX_MODEL_H

#include <cmath>
#include <cassert>

namespace prox
{

  namespace detail
  {

    /**
     *
     * This function computes the proximal point of a 2D box, meant to emulate the NCP friction formulation.
     */
    template <typename T>
    inline static void box_model(
                                 T const & z_s
                                 , T const & z_t
                                 , T const & z_tau
                                 , T const & mu_s
                                 , T const & mu_t
                                 , T const & mu_tau
                                 , T const & lambda_n
                                 , T & lambda_s
                                 , T & lambda_t
                                 , T & lambda_tau
                                 )
    {

      using std::min;
      using std::max;

      T const n_s   = mu_s* max<T>( lambda_n, 0 );
      T const n_t   = mu_t* max<T>( lambda_n, 0 );

      assert( is_number( n_s )    || !"box_model(): a was not a number");
      assert( n_s >= 0    || !"box_model(): a non-positive");
      assert( is_number( n_t )    || !"box_model(): a was not a number");
      assert( n_t >= 0    || !"box_model(): a non-positive");

      lambda_s      =    min( n_s, max(z_s, -n_s ) );
      lambda_t      =    min( n_t, max(z_t, -n_t ) );
      lambda_tau    =    0;

    }

  } // namespace detail
} // namespace prox

// PROX_BOX_MODEL_H
#endif
