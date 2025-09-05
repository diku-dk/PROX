#ifndef PROX_GET_EXTERNAL_FORCES_VECTOR_H
#define PROX_GET_EXTERNAL_FORCES_VECTOR_H

#include "prox_math_policy.h"
#include "prox_update_inertia_tensor.h"
#include <prox_force_callbacks.h>

#include <tiny_is_number.h>
#include <tiny_is_finite.h>

#include <cassert>

namespace prox
{

  /**
   *
   * @param h      Upon return this parameter contains the total external forces
   *               and torques acting on the bodies in the system
   */
template <typename body_iterator, typename T>
inline void get_external_forces_vector(body_iterator begin, body_iterator end,
                                       Gravity< T > const& gravity, Damping< T > const& damping,
                                       typename prox::MathPolicy<T>::vector6_type& h)
{
    typedef typename prox::MathPolicy<T>::block6x1_type B6x1;

    size_t const N = std::distance(begin,end);
    h.resize( N );

    size_t k = 0u;
    for(body_iterator body = begin;body!=end;++body, ++k)
    {
      B6x1 & b = h( k );

      if(body->is_fixed() || body->is_scripted() )
      {
        b(0) = 0;
        b(1) = 0;
        b(2) = 0;
        b(3) = 0;
        b(4) = 0;
        b(5) = 0;
        continue;
      }

      EigenVector3<T> total_force(0, 0, 0);
      auto total_torque = total_force;
      auto force = total_force;
      auto torque = total_torque;

      //--- First we add global world gravity force ----------------------------
      gravity.compute_force_and_torque((*body), force, torque);

      total_force  += force;
      total_torque += torque;

      //--- Second we add global world damping force ---------------------------
      damping.compute_force_and_torque((*body), force, torque);

      total_force  += force;
      total_torque += torque;

      //--- Third we add any local body forces that might be applied -----------
      for (const auto& callback : body->get_force_callbacks())
      {
          callback->compute_force_and_torque((*body), force, torque);
          total_force += force;
          total_torque += torque;
      }

      assert(is_number(total_force(0)) || !"get_external_forces_vector(): Nan");
      assert(is_number(total_force(1)) || !"get_external_forces_vector(): Nan");
      assert(is_number(total_force(2)) || !"get_external_forces_vector(): Nan");

      assert(is_number(total_torque(0)) || !"get_external_forces_vector(): Nan");
      assert(is_number(total_torque(1)) || !"get_external_forces_vector(): Nan");
      assert(is_number(total_torque(2)) || !"get_external_forces_vector(): Nan");

      assert(is_finite(total_force(0)) || !"get_external_forces_vector(): Inf");
      assert(is_finite(total_force(1)) || !"get_external_forces_vector(): Inf");
      assert(is_finite(total_force(2)) || !"get_external_forces_vector(): Inf");

      assert(is_finite(total_torque(0)) || !"get_external_forces_vector(): Inf");
      assert(is_finite(total_torque(1)) || !"get_external_forces_vector(): Inf");
      assert(is_finite(total_torque(2)) || !"get_external_forces_vector(): Inf");

      b(0) = total_force(0);
      b(1) = total_force(1);
      b(2) = total_force(2);

      auto const& w = body->get_spin();
      auto const& I_bf = body->get_inertia_bf();

      assert(is_number(w(0))|| !"get_external_forces_vector(): Nan");
      assert(is_number(w(1))|| !"get_external_forces_vector(): Nan");
      assert(is_number(w(2))|| !"get_external_forces_vector(): Nan");

      assert(is_finite(w(0))|| !"get_external_forces_vector(): Inf");
      assert(is_finite(w(1))|| !"get_external_forces_vector(): Inf");
      assert(is_finite(w(2))|| !"get_external_forces_vector(): Inf");

      auto R = body->get_orientation().toRotationMatrix();

      EigenMatrix3<T> I{
          {1, 0, 0},
          {0, 1, 0},
          {0, 0, 1}
      };
      detail::update_inertia_tensor(R, I_bf, I);

      auto const wIw = w.cross(I * w);

      assert(is_number(wIw(0)) || !"get_external_forces_vector(): Nan");
      assert(is_number(wIw(1)) || !"get_external_forces_vector(): Nan");
      assert(is_number(wIw(2)) || !"get_external_forces_vector(): Nan");

      assert(is_finite(wIw(0)) || !"get_external_forces_vector(): Inf");
      assert(is_finite(wIw(1)) || !"get_external_forces_vector(): Inf");
      assert(is_finite(wIw(2)) || !"get_external_forces_vector(): Inf");

      b(3) = total_torque(0) - wIw(0);
      b(4) = total_torque(1) - wIw(1);
      b(5) = total_torque(2) - wIw(2);
    }
}
}// namespace prox

// PROX_GET_EXTERNAL_FORCES_VECTOR_H
#endif
