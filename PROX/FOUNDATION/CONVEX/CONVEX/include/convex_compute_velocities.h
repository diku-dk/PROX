#ifndef CONVEX_COMPUTE_VELOCITIES_H
#define CONVEX_COMPUTE_VELOCITIES_H

#include <cmath>
#include <cassert>

#include <tiny_math_types.h>

namespace convex
{

  /**
   * Compute Velocities.
   * This function tries to convert two poses from the motion of an object into equivalent velocities.
   *
   * The function works under the assumption that the object in questions moves at constant
   * linear and angular velocities between the two given poses. The duration of the motion
   * is assumed to be one unit (second).
   *
   * @param T_from    A coordinate transformation indicating the starting pose of the motion.
   * @param T_to      A coordinate transformation indicating the ending pose of the motion.
   * @param delta_tau The time between from pose and to pose.
   * @param v         Upon return this argument holds the value of the constant linear velocity.
   * @param omega     Upon return this argument holds the value of the constant angular velocity.
   */
  template< typename T>
  inline void compute_velocities(const CoordSysEigen<T>& T_from
                          , const CoordSysEigen<T>& T_to
                          , const T& delta_tau
                          , EigenVector3<T>& v
                          , EigenVector3<T> & omega
                          )
  {
    using std::atan2;


    assert(  delta_tau > 0 || !"compute_velocities(): time step must be positive");

    // Translation is straightforward
    v = ((T_to.T() - T_from.T()) / delta_tau);

    T theta;
    EigenVector3<T> n;
    getAxisAngle(
        ( T_to.Q()*( T_from.Q() ).conjugate() )   // Change in orientation from ``from'' to ''to'', ie. R = T_to * T_from^{-1}
                         , n
                         , theta
                         );
    omega = ((theta/ delta_tau)*n);
  }

} // namespace convex

// CONVEX_COMPUTE_VELOCITIES_H
#endif
