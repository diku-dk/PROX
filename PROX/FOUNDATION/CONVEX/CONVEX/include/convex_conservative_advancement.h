#ifndef CONVEX_CONSERVATIVE_ADVANCEMENT_H
#define CONVEX_CONSERVATIVE_ADVANCEMENT_H

#include <cassert>
#include <eigenhelperall.h>
#include <convex_compute_closest_points.h>
#include <convex_integrate_motion.h>

namespace convex
{

  /**
   * Conservative Advancement.
   * This function tries to determine whether two objects have impacted
   * during their motion. The function assumes that the objects will move
   * continuously in the future with contact linear and angular velocities.
   *
   *
   * @param r_A          Initial position of shape A.
   * @param q_A          Initial quaternion orientation of shape A.
   * @param v_A          Linear velocity of shape A.
   * @param w_A          Angular velocity of shape A.
   * @param A            Pointer to convex shape A.
   * @param r_max_A      Maximum radius of the shape of object A.
   *
   * @param r_B          Initial position of shape B.
   * @param q_B          Initial quaternion orientation  of shape B.
   * @param v_B          Linear velocity of shape B.
   * @param w_B          Angular velocity of shape B.
   * @param B            Pointer to convex shape B.
   * @param r_max_B      Maximum radius of the shape of object B.
   *
   * @param p_A             Upon return this argument holds the cloest point on object A in case of an impact.
   * @param p_B             Upon return this argument holds the cloest point on object B in case of an impact.
   * @param time_of_impact  Upon return if an impact is found then this argument holds the estimated value of the time of impact.
   * @param iterations      Upon return this argument holds the number of used iterations by the function. If the value is equal to the max_iterations argument then the function did not converge to an answer.
   * @param epsilon         The size of the collision envelope. That is the smallest separation distance between A and B where we consider A and B to be in touching contact.
   * @param max_tau         The maximum time into the future that the function will look for a time of impact.
   * @param max_iterations  The maximum number of allowed iterations that the function can take.
   *
   * @return                If an impact is found then the return value is true otherwise it is false.
   */
  template<typename T>
  inline bool conservative_advancement(const CoordSysEigen<T>& X_A
                                , const EigenVector3<T>& v_A
                                , const EigenVector3<T>& w_A
                                , geometry::SupportMapping<T> const * A
                                , const T r_max_A
                                , const CoordSysEigen<T>& X_B
                                , const EigenVector3<T>& v_B
                                , const EigenVector3<T>& w_B
                                , geometry::SupportMapping<T> const * B
                                , const T& r_max_B
                                , EigenVector3<T>& p_A
                                , EigenVector3<T>& p_B
                                , T& time_of_impact
                                , size_t & iterations
                                , const T& epsilon
                                , const T& max_tau
                                , size_t const & max_iterations
                                )
  {

    assert( r_max_A > 0              || !"conservative_advancement(): maximum distance of object A must be positive");
    assert( r_max_B > 0              || !"conservative_advancement(): maximum distance of object B must be positive");
    assert( max_tau > 0              || !"conservative_advancement(): maximum time-step must be positive");
    assert( epsilon > 0              || !"conservative_advancement(): collision envelope must be positive");
    assert( max_iterations > 0u               || !"conservative_advancement(): maximum iterations must be positive");
    assert(epsilon >= (1e-2)  || !"conservative_advancement(): Too aggressive setting of epsilon, compute_closest_points uses tolerance 10e4");

    T tau = 0;

    for(iterations=1u; iterations <= max_iterations; ++iterations)
    {
      // Compute the coordinate transformations corresponding to the current tau value
        CoordSysEigen<T> T_A = integrate_motion<T>( X_A, tau, (v_A), (w_A) );
        CoordSysEigen<T> T_B = integrate_motion<T>( X_B, tau, (v_B), (w_B) );

      // Compute the closest points at the time tau
      compute_closest_points<T>( (T_A), A, (T_B), B, (p_A), (p_B) );

      // Estimate normal direction and current minimum distance between A and B
      EigenVector3<T> v = (p_A) - (p_B);

      // 2015-11-19 Kenny: If GJK did not converge completely then p_A and p_B will
      //                   be slightly off... this means the distance we compute here
      //                   is in fact a little larger than the "true" minimum distance. We
      //                   counter this by using the treshold test below. However, if one
      //                   is too aggressive then one might overstep the true TOI estimate
      //                   by a tiny fraction.

      T min_distance = (v).norm();

      if( min_distance <= epsilon )
      {
        time_of_impact = tau;
        return true;
      }

      // 2015-11-19 Kenny: If GJK did not converge completely then p_A and p_B will
      //                   be slightly off... this means the "n" direction can be off. The
      //                   problem is worsen as objects come close, as n will be determined
      //                   by "substracting" numbers that gets smaller and smaller....
      //
      //                   At some point the impression will dominate and determine the
      //                   direction more than the "closes" point will dominate the
      //                   direction.
      //
      //                   One should use epislon to avoid getting into this case... too
      //                   agressive epsilon means that max_velocity will be badly estimated
      //                   (too low than the true value.
      //

      EigenVector3<T> n = ( v ).normalized();

      // Estimate maximum relative normal velocity between any two points from A and B

      T max_velocity = (v_B - v_A).dot(n) + (w_A).norm()*r_max_A + (w_B).norm()*r_max_B;

      if (max_velocity <= 0 )
        return false;

      // Compute conservative lower bound for when A and B could impact
      T delta_tau = min_distance / max_velocity;
      tau += delta_tau;

      if ( tau > max_tau )
        return false;
    }

    // not enough iterations to determine what goes on! We give up
    return false;
  }

} // namespace convex

// CONVEX_CONSERVATIVE_ADVANCEMENT_H
#endif
