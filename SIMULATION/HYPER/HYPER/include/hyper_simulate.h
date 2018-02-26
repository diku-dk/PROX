#ifndef HYPER_SIMULATE_H
#define HYPER_SIMULATE_H

#include <hyper_engine.h>
#include <hyper_collision_detection.h>
#include <hyper_compute_gravity_forces.h>
#include <hyper_semi_implicit_time_step.h>
#include <hyper_adaptive_time_step.h>
#include <hyper_modifiers.h>
#include <hyper_compute_CFL_time_step_size.h>

#include <cmath>
#include <cassert>

namespace hyper
{

  /**
   * Simulate Engine Forward in Time.
   *
   * @param engine   The engine holding the configuration that should be time-integrated forward in time.
   * @param time     The amount of time to simulate forward in time.
   */
  template<typename MT>
  inline void simulate( Engine<MT> & engine, typename MT::real_type const & time )
  {
    using std::min;

    typedef typename MT::real_type       T;
    typedef typename MT::vector3_type    V;
    typedef typename MT::value_traits    VT;
    typedef typename MT::quaternion_type Q;

    typedef          Engine<MT>                            engine_type;
    typedef typename engine_type::body_type                body_type;
    typedef typename engine_type::body_iterator            body_iterator;
    typedef typename engine_type::const_body_iterator      const_body_iterator;
    typedef typename engine_type::contact_type             contact_type;
    typedef typename engine_type::contact_iterator         contact_iterator;
    typedef typename engine_type::const_contact_iterator   const_contact_iterator;
    typedef typename engine_type::params_type              params_type;


    assert( time > VT::zero()                        || !"simulate(): time to simulate forward must be positive");
    assert( engine.params().time_step() > VT::zero() || !"simulate(): time to simulate forward must be positive");

    if( norm(engine.gravity())> VT::zero() )
    {
      for (body_iterator body = engine.body_begin();body != engine.body_end(); ++body)
      {
        compute_gravity_forces(
                               body->m_mesh
                               , body->m_neighbors
                               , engine.gravity()
                               , body->m_x0
                               , body->m_model
                               , body->m_Fext
                               );

      }
    }

    T const time_step = min( engine.params().time_step(), time );

    T time_left = time;

    while(time_left > VT::zero())
    {
      T const dt_wanted  = min(time_left, time_step);
      T const dt_CFL     = engine.params().use_cfl_condition() ? compute_CFL_time_step_size( engine, dt_wanted ) : dt_wanted;
      T const dt         = dt_CFL;

      for (body_iterator body = engine.body_begin();body != engine.body_end(); ++body)
        MT::convert(body->m_x, body->m_X, body->m_Y, body->m_Z);

      collision_detection( engine );

      for (body_iterator body = engine.body_begin();body != engine.body_end(); ++body)
      {

        if(body->m_scripted_motion)
        {
          body->m_scripted_motion->compute( *body, engine.time() );
        }
        else
        {
          compute_traction_forces(body->m_traction_conditions, body->m_mesh, body->m_x, body->m_F);

          sparse::add(body->m_Fext, body->m_F);

          switch (engine.params().get_time_step_method_type())
          {
            case params_type::semi_implicit_type:

              semi_implicit_time_step(
                                      body->m_mesh
                                      , body->m_neighbors
                                      , body->m_dirichlet_conditions
                                      , dt
                                      , body->m_x
                                      , body->m_x0
                                      , body->m_v
                                      , body->m_F
                                      , body->m_model
                                      , jacobi_precond
                                      );

              break;

            case params_type::adaptive_type:

              adaptive_time_step(dt, *body, engine.params() );

              break;

              
            default:
              break;
          }



        }

        MT::convert(body->m_x, body->m_X, body->m_Y, body->m_Z);
      }

      time_left     -= dt;
      engine.time() += dt;
    }

  }

} // namespace hyper

// HYPER_SIMULATE_H
#endif
