#ifndef PROX_POSITION_UPDATE_H
#define PROX_POSITION_UPDATE_H

#include "prox_math.h"
#include "tiny_math_types.h"
#include <cmath>
#include <stdexcept>

namespace prox
{

template <typename T>
void position_update(const NCVec7<T>& q, const NCVec6<T>& u, T dt, NCVec7<T>& qnew)
{
    size_t const N = u.size();

    if( q.size() != N )
      throw std::logic_error("position_update(): s has incorrect dimension");

    // 2009-08-04 Kenny: You are right this seems pointless. In fact I think the correct behaviour should be to get the number
    // of blocks from q, then verify if u has the same number of blocks, and lastly if qnew is big enough to hold the results
    if( u.size() != N )// 2009-08-03 Sarah: is this not redundant since N = u.size()?
      throw std::logic_error("position_update(): u has incorrect dimension");

    if(&q != &qnew)
    {
      qnew.resize( N );
    }

    float const dt_half = dt / 2;

    for(size_t i = 0; i<N; ++i)
    {
        auto& qnew_b = qnew(i);
        auto const& q_b = q(i);
        auto const& u_b = u(i);

        qnew_b(0) = q_b(0) + dt * u_b(0);
        qnew_b(1) = q_b(1) + dt * u_b(1);
        qnew_b(2) = q_b(2) + dt * u_b(2);

        EigenQuaternion<T> Q;
        Q.w() = q_b(3);
        Q.x() = q_b(4);
        Q.y() = q_b(5);
        Q.z() = q_b(6);

        EigenVector3<T> W;
        W(0) = u_b(3);
        W(1) = u_b(4);
        W(2) = u_b(5);

        bool const finite = true; // Kenny: Hmm, this might have to be a user-specified paramter?

        if (finite)
        {
            // Do a finitedimensional update instead
            //
            // Create the quaternion, H, corresonding to a rotation of
            //
            //    \theta = \Delta t \norm{\vec \omega}
            //
            // Around the unit-axis
            //
            //     \vec n = \frac{\vec \omega}{  \norm{\vec \omega} }
            //
            // That is
            //
            //    H = [ \cos\left( \frac{\theta}{2} \right), \sin\left( \frac{\theta}{2} \right) \vec n ]
            //
            // Now perform the update of the orientation simply by
            //
            //   Q \leftarrow  H Q
            //
            // There should be no need to normalize Q after this operation.
            float const radian = norm(W) * dt;
            auto axis = W.normalized();
            auto R = Rotateu(radian, axis);
            Q = R * Q;
        }
      else
      {
          //--- Just do an infinitesimal rotation update (i.e. a forward euler step)
          // Compute W*Q
          EigenVector3<T> qimag(Q.x(), Q.y(), Q.z());
          auto newReal = -dot(W, qimag);
          EigenVector3<T> newImag = cross(W, qimag) + W * Q.w();
          Q.w() += newReal * dt * .5f;
          Q.x() += newImag.x() * dt * .5f;
          Q.y() += newImag.y() * dt * .5f;
          Q.z() += newImag.z() * dt * .5f;
      }

      //--- To counter-act numerical problems
      Q = Q.normalized();

      qnew_b(3) = Q.w();
      qnew_b(4) = Q.x();
      qnew_b(5) = Q.y();
      qnew_b(6) = Q.z();
    }
}
} // namespace prox

// PROX_POSITION_UPDATE_H
#endif
