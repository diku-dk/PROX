#ifndef PROX_STEPPER_BIND_H
#define PROX_STEPPER_BIND_H

#include <prox_enums.h>
#include <prox_stepper.h>
#include <prox_stepper_moreau.h>
#include <prox_stepper_semi_implicit.h>
#include <prox_stepper_empty.h>

#include <broad.h>
#include <narrow.h>

#include <util_log.h>

#include <cassert>
#include <vector>

namespace prox
{
  
  /**
   * Solver function pointer binder.
   */
  template< typename M  >
  class StepperBinder
  : public Stepper< M >
  {
  public:
    
    typedef void func_type(
                           typename M::real_type const & dt
                           , std::vector< RigidBody< M > > & bodies
                           , std::vector< std::vector< Property< M > > > const &  properties
                           , Gravity< M > const & gravity
                           , Damping< M > const & damping
                           , Params< M > const & params
                           , broad::System<typename M::real_type> & broad_system
                           , narrow::System<typename M::tiny_types> & narrow_system
                           , std::vector< ContactPoint< M > > & contacts
                           , M const & tag
                           );
    
    func_type * m_stepper;
    
  public:
    
    StepperBinder()
    : m_stepper(0)
    {
    }
    
    StepperBinder(func_type * stepper)
    : m_stepper(stepper)
    {
    }
    
  public:

    void operator()(
                      typename M::real_type const & dt
                    , std::vector< RigidBody< M > > & bodies
                    , std::vector< std::vector< Property< M > > > const & properties
                    , Gravity< M > const & gravity
                    , Damping< M > const & damping
                    , Params< M > const & params
                    , broad::System<typename M::real_type> & broad_system
                    , narrow::System<typename M::tiny_types> & narrow_system
                    , std::vector< ContactPoint< M > > & contacts
                    , M const & tag
                    ) const
    {
      assert( this->m_stepper || !"StepperBinder(): stepper was null");
      
      this->m_stepper( dt, bodies, properties, gravity, damping, params, broad_system, narrow_system, contacts, tag );
    }
    
  };
  
  /**
   *
   */
  template< typename M  >
  inline StepperBinder<M> bind_stepper( stepper_type const & type )
  {
    util::Log logging;

    switch( type )
    {
      case moreau:
        logging << "bind_stepper(): using moreau stepper"<< util::Log::newline();
        return StepperBinder<M>( &moreau_time_stepper<M> );
        
      case semi_implicit:
        logging << "bind_stepper(): using semi-implicit stepper"<< util::Log::newline();
        return StepperBinder<M>( &semi_implicit_time_stepper<M> );

      case empty:
        logging << "bind_stepper(): using empty stepper"<< util::Log::newline();
        return StepperBinder<M>( &empty_stepper<M> );
        
      default:
        assert(!"bind_stepper(): unknown stepper type");
        break;
    };
    
    return StepperBinder<M>();
  }
  
} //namespace prox

// PROX_STEPPER_BIND_H
#endif