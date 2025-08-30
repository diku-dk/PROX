#ifndef PROX_ANALYTICAL_ELLIPSOID_H
#define PROX_ANALYTICAL_ELLIPSOID_H

#include <tiny_is_number.h>
#include <tiny_polynomial_roots.h>
#include <tiny.h>

namespace prox
{
  namespace detail
  {
    /**
    * This function is for isotropic planar friction.
    *
    * This function sets up a polynomial and solves analytically
    * for the roots of the polynomial.
    */
    template <typename T>
    inline static void analytical_ellipsoid(
                                                 T const & z_s
                                                 , T const & z_t
                                                 , T const & z_tau
                                                 , T const & mu
                                                 , T const & mu_tau
                                                 , T const & lambda_n
                                                 , T & lambda_s
                                                 , T & lambda_t
                                                 , T & lambda_tau

      )
    {
      typedef tiny::ValueTraits<T>   value_traits;

      if(lambda_n <= 0)
      {
        lambda_s   = 0;
        lambda_t   = 0;
        lambda_tau = 0;
        return;
      }

      T const a = mu*lambda_n;
      T const b = a;
      T const c = mu_tau*lambda_n;

      assert( is_number( a ) || !"analytical_ellipsoid(): a was not a number");
      assert( is_number( b ) || !"analytical_ellipsoid(): b was not a number");
      assert( is_number( c ) || !"analytical_ellipsoid(): c was not a number");
      assert( a > 0       || !"analytical_ellipsoid(): a non-positive");
      assert( b > 0       || !"analytical_ellipsoid(): a non-positive");
      assert( c > 0       || !"analytical_ellipsoid(): a non-positive");

      // let the ellipsoid surface be given by the zero-level set contour of the function, f(z) = 0
      //
      //   f(z) = z_s^2/a^2 + z_t^2/b^2 + z_tau^2/c^2 - 1 = 0
      //
      //
      T const e = ((z_s*z_s) / (a*a)) + ((z_t*z_t) / (b*b)) + ((z_tau*z_tau) / (c*c));

      if( e <= 1 )
      {
        lambda_s   = z_s;
        lambda_t   = z_t;
        lambda_tau = z_tau;
        return;
      }

      // To find the closest point, lambda, on the surface of the ellipsoid
      // we are seeking the value of the line paramter, alpha, such that
      //
      //    lambda + alpha nabla f(lambda) = z  (*1)
      //
      // where alpha >= 0 and the outward gradient is given by
      //
      //                  | 2 z_s /a^2   |
      //    nabla f(z) =  | 2 z_t /b^2   |
      //                  | 2 z_tau /c^2 |
      //
      // Thus from (*1) we have
      //
      //   lambda = (I + 2 alpha diag( 1/a^2, 1/b^2, 1/b^2 ))^{-1} z
      //
      // using this to elimate lambda in f(lambda) = 0 and
      // assuming b= a yields fourth order polynomial
      //
      //    lambda_s   =   a^2 z_s / (a^2 + 2 alpha)
      //    lambda_t   =   a^2 z_t / (a^2 + 2 alpha)
      //    lambda_tau =  c^2 z_tau / (c^2 + 2 alpha)
      //
      //   f(lambda) = a^2 (z_s^2+z_t^2) / (a^2 + 2 alpha)^2 + c^2 (z_tau^2) / (c^2 + 2 alpha)^2 - 1= 0
      //
      // making the subsitution t = 2 alpha we have
      //
      //
      //   (a^2 + t)^2(c^2 + t)^2  - a^2 (z_s^2+z_t^2)(a^2 + t)^2 -  c^2 (z_tau^2)(a^2 + t)^2 = 0
      //
      //      a4 t^4 + a3 t^3 + a2 t^2 + a1 t^1 + a0 = 0
      // where
      //
      //    a4 = 1
      //    a3 = 2(c^2 + a^2)
      //    a2 = c^4 + a^4 + 4 a^2 c^2 - c^2 z_tau^2 - a^2 (z_s^2 + z_t^2)
      //    a1 = 2( a^4 c^2 + c^4 a^2 - c^2 a^2 ( z^T z ))
      //    a0 = a^4 c^4 - c^2 a^4 z_tau^2 - a^2 c^4 (z_s^2 + z_t^2)
      //
      // See also PhD thesis by Studer on page 182

      T const a2 = a*a;
      T const c2 = c*c;
      T const a4 = a2*a2;
      T const c4 = c2*c2;
      T const x2 = z_s*z_s;
      T const y2 = z_t*z_t;
      T const z2 = z_tau*z_tau;

      T const coef4 = 1;
      T const coef3 = 2*(c2+a2);
      T const coef2 = c4 + a4 + 4*a2*c2 - c2*z2 - a2*(x2 + y2);
      T const coef1 = 2*(a4*c2 + c4*a2 - c2*a2*(x2+y2+z2));
      T const coef0 = a4*c4 - c2*a4*z2 - a2*c4*(x2+y2);

      // We will have two imaginary roots, one negative and one positive
      // root. Only the positive root is of interest.
      T t = 0;
      {
        unsigned int count = 0u;

        T roots[4] = {0, 0, 0, 0};

        tiny::compute_polynomial_roots(coef0,coef1,coef2,coef3,coef4,count, roots );

        for(size_t i =0;i < count;++i)
          t = ( roots[i] > t ) ? roots[i] : t;
      }

      // Knowing the alpha-value we can compute the closest point using (*1)
      lambda_s =  z_s / (1 + t/(a*a));
      lambda_t =  z_t / (1 + t/(b*b));
      lambda_tau =  z_tau / (1 + t/(c*c));
    }


    /**
     *
     * Overloaded version. The only purpose of this version is
     * such that all sub solvers have the same function signature.
     */
    template <typename T>
    inline static void analytical_ellipsoid(
                                              T const & z_s
                                              , T const & z_t
                                              , T const & z_tau
                                              , T const & mu_s
                                              , T const & /*mu_t*/
                                              , T const & mu_tau
                                              , T const & lambda_n
                                              , T & lambda_s
                                              , T & lambda_t
                                              , T & lambda_tau
                                              )
    {
      analytical_ellipsoid(z_s, z_t, z_tau, mu_s, mu_tau,lambda_n,lambda_s,lambda_t, lambda_tau);
    }


  } // namespace detail
} // namespace prox

// PROX_ANALYTICAL_ELLIPSOID_H
#endif
