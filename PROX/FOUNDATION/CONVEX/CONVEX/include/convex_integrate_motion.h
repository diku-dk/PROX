#ifndef CONVEX_INTEGRATE_MOTION_H
#define CONVEX_INTEGRATE_MOTION_H

#include <cmath>
#include <cassert>

#include <eigenhelperall.h>

namespace convex
{

  /**
   * Integrate Motion.
   *
   * @param X      The current coordinate transformation of the object.
   * @param tau    The time step into the future where the coordinate transformation of the object should be computed.
   * @param v      The current linear velocity of the object.
   * @param omega  The current angular velocity of the object.
   *
   * @return       Upon return this argument holds the coordinate transformation of the object at time tau.
   */
  template< typename T>
  inline CoordSysEigen<T> integrate_motion(
                                             const CoordSysEigen<T>& X
                                             , const T tau
                                             , const EigenVector3<T>& v
                                             , const EigenVector3<T>& omega
                                             )
  {


    assert( tau >= 0 || !"integrate_motion(): Tau must be non-negative");

      T const radian           = tau * ( omega ).norm();
    const EigenVector3<T> axis             = ( omega ).normalized();

    assert( is_number( radian )  || !"integrate_motion(): NaN encountered");
    assert( is_number( axis(0) ) || !"integrate_motion(): NaN encountered");
    assert( is_number( axis(1) ) || !"integrate_motion(): NaN encountered");
    assert( is_number( axis(2) ) || !"integrate_motion(): NaN encountered");

    EigenQuaternion<T> dq;
    EigenVector3<T> dv;
    dq = Rotateu( radian, axis);
    dv = v*tau;

    assert( is_number( dv(0) ) || !"integrate_motion(): NaN encountered");
    assert( is_number( dv(1) ) || !"integrate_motion(): NaN encountered");
    assert( is_number( dv(2) ) || !"integrate_motion(): NaN encountered");

    return CoordSysEigen<T>( dv + X.T(), (dq*X.Q()) );
  }

} // namespace convex

// CONVEX_INTEGRATE_MOTION_H
#endif
