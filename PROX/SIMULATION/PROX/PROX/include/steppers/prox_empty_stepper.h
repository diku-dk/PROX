#ifndef PROX_EMPTY_TIME_STEPPER_H
#define PROX_EMPTY_TIME_STEPPER_H

#include <prox_rigid_body.h>
#include <prox_contact_point.h>

#include <prox_collision_detection.h>

#include <prox_params.h>
#include <prox_math_policy.h>

#include <broad.h>
#include <narrow.h>

#include <util_profiling.h>
#include <util_log.h>

namespace prox
{
  /**
   * Empty time stepper.
   * This essentially does nothing, except for invoking the collision detection system.
   */
template < typename T>
inline void empty_stepper(T dt, std::vector< RigidBody<T> >& bodies,
                          std::vector< std::vector< Property<T> > > const& properties,
                          Gravity< T > const& gravity, Damping< T > const& damping,
                          Params<T> const& params, broad::System<T>& broad_system,
                          narrow::System<T>& narrow_system,
                          std::vector< ContactPoint<T> >& contacts)
{
    using M = tiny::MathTypes<T>;
    M tag;
    util::Log logging;

    START_TIMER("stepper");

    collision_detection(
                        bodies
                        , broad_system
                        , narrow_system
                        , contacts
                        , params
                        , tag
                        );

    unsigned int const number_of_contacts = contacts.size();

    logging << "empty_stepper(): Number of contacts = " << number_of_contacts << util::Log::newline();

    STOP_TIMER("stepper");
}

} //namespace prox

// PROX_EMPTY_TIME_STEPPER_H
#endif
