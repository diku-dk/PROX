#ifndef PROX_GJK_ELLIPSOID_H
#define PROX_GJK_ELLIPSOID_H

#include <geometry.h>   // needed for geometry::Point
#include <convex.h>
#include <cassert>

namespace prox
{
    namespace detail
    {

        /**
         *
         * Complete general purpose force-law.
         *
         *
         * This implementation uses the CONVEX algorithm for
         * determining the proximal point.
         *
         * This is actually pretty cool, because if we in the
         * future want non-ellipsoid contact laws then the
         * implementation already supports this!
         */
        template <typename T>
        inline static void gjk_ellipsoid(
                                              T const & z_s,
                                              T const & z_t,
                                              T const & z_tau,
                                              T const & mu_s,
                                              T const & mu_t,
                                              T const & mu_tau,
                                              T const & lambda_n,
                                              T & lambda_s,
                                              T & lambda_t,
                                              T & lambda_tau
                                              )
        {

            if ( lambda_n <= 0 )
            {
                lambda_s   = 0;
                lambda_t   = 0;
                lambda_tau = 0;
                return;
            }

            T const a = mu_s*lambda_n;
            T const b = mu_t*lambda_n;
            T const c = mu_tau*lambda_n;

            assert( is_number( a ) || !"gjk_ellipsoid(): a was not a number");
            assert( is_number( b ) || !"gjk_ellipsoid(): b was not a number");
            assert( is_number( c ) || !"gjk_ellipsoid(): c was not a number");
            assert( a > 0       || !"gjk_ellipsoid(): a non-positive");
            assert( b > 0       || !"gjk_ellipsoid(): a non-positive");
            assert( c > 0       || !"gjk_ellipsoid(): a non-positive");

            geometry::Point<T>                       point;
            convex::Ellipsoid<T>  ellipsoid;

            point.coord() = EigenVector3<T>( z_s, z_t, z_tau);

            // TODO: Chek theory
            ellipsoid.setScale({a,b,c});  // TODO check this is how to setup the scale!

            size_t    const max_iterations       = 100u;
            T const absolute_tolerance = T(10e-6);
            T const relative_tolerance = T(10e-6);
            T const stagnation_tolerance = T(10e-15);

            CoordSysEigen<T> transformA;
            CoordSysEigen<T> transformB;
            transformA.T() = {0,0,0};
            transformA.Q() = EigenQuaternion<T>::Identity();//may cause issues on Windows/VS platform
            transformB.T() = {0,0,0};
            transformB.Q() = EigenQuaternion<T>::Identity();

            EigenVector3<T> pa;
            EigenVector3<T> pb;
            size_t iterations     = 0u;
            size_t status         = 0u;
            T distance    = std::numeric_limits<T>::max();

            convex::compute_closest_points<T>(
                                                                            transformA
                                                                            , &point
                                                                            , transformB
                                                                            , &ellipsoid
                                                                            , pa
                                                                            , pb
                                                                            , distance
                                                                            , iterations
                                                                            , status
                                                                            , absolute_tolerance
                                                                            , relative_tolerance
                                                                            , stagnation_tolerance
                                                                            , max_iterations
                                                                            );

            assert( status != convex::ITERATING                     || !"gjk_ellipsoid(): gjk internal error");
            assert( status != convex::EXCEEDED_MAX_ITERATIONS_LIMIT || !"gjk_ellipsoid(): gjk internal error");
            assert( status != convex::NON_DESCEND_DIRECTION         || !"gjk_ellipsoid(): gjk internal error");

            lambda_s = pb(0);
            lambda_t = pb(1);
            lambda_tau = pb(2);
        }

    } // namespace detail
} // namespace prox

// PROX_GJK_ELLIPSOID_H
#endif
