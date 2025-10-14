#ifndef PROX_POSITION_UPDATE_H
#define PROX_POSITION_UPDATE_H

#include "eigenhelperfunctions.h"
#include <cmath>
#include <stdexcept>

namespace prox
{

template <typename T>
void position_update_eigen(const Eigen::VectorX<T>& q,
                           const Eigen::VectorX<T>& u, T dt,
                           Eigen::VectorX<T>& qnew)
{
    size_t const N = u.size() / 6; // Each body has 6 DOF in u

    if (q.size() != N * 7)
    {
        throw std::logic_error("position_update(): q has incorrect dimension");
    }

    if (u.size() != N * 6)
    {
        throw std::logic_error("position_update(): u has incorrect dimension");
    }

    if (qnew.size() != N * 7) { qnew.resize(N * 7); }

    T const dt_half = dt / 2;

    for (size_t i = 0; i < N; ++i)
    {
        // Calculate indices for this body
        size_t q_idx = i * 7;
        size_t u_idx = i * 6;
        size_t qnew_idx = i * 7;

        // Update position
        qnew(qnew_idx) = q(q_idx) + dt * u(u_idx);
        qnew(qnew_idx + 1) = q(q_idx + 1) + dt * u(u_idx + 1);
        qnew(qnew_idx + 2) = q(q_idx + 2) + dt * u(u_idx + 2);

        // Extract quaternion from q
        Eigen::Quaternion<T> Q;
        Q.w() = q(q_idx + 3);
        Q.x() = q(q_idx + 4);
        Q.y() = q(q_idx + 5);
        Q.z() = q(q_idx + 6);

        // Extract angular velocity from u
        Eigen::Matrix<T, 3, 1> W;
        W(0) = u(u_idx + 3);
        W(1) = u(u_idx + 4);
        W(2) = u(u_idx + 5);

        bool const finite = true;

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
            T const radian = W.norm() * dt;
            if (radian > std::numeric_limits<T>::epsilon())
            {
                Eigen::Matrix<T, 3, 1> axis = W.normalized();
                Eigen::AngleAxis<T> R(radian, axis);
                Q = R * Q;
            }
        }
        else
        {
            // Infinitesimal rotation update
            Eigen::Matrix<T, 3, 1> qimag(Q.x(), Q.y(), Q.z());
            T newReal = -W.dot(qimag);
            Eigen::Matrix<T, 3, 1> newImag = W.cross(qimag) + W * Q.w();

            Q.w() += newReal * dt_half;
            Q.x() += newImag.x() * dt_half;
            Q.y() += newImag.y() * dt_half;
            Q.z() += newImag.z() * dt_half;
        }

        // Normalize quaternion
        Q.normalize();

        // Store updated quaternion
        qnew(qnew_idx + 3) = Q.w();
        qnew(qnew_idx + 4) = Q.x();
        qnew(qnew_idx + 5) = Q.y();
        qnew(qnew_idx + 6) = Q.z();
    }
}
} // namespace prox

// PROX_POSITION_UPDATE_H
#endif
