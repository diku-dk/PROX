#ifndef PROX_BIND_STEPPER_H
#define PROX_BIND_STEPPER_H

#include <prox_enums.h>

#include <steppers/prox_stepper.h>

#include <steppers/prox_moreau_time_stepper.h>
#include <steppers/prox_semi_implicit_time_stepper.h>
#include <steppers/prox_empty_stepper.h>

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
template <typename T> class StepperBinder
{
public:
    StepperBinder(stepper_type type = moreau)
        : m_type{type}
    {
        util::Log logging;
        switch (type)
        {
        case moreau :        logging << "bind_stepper(): using moreau stepper" << util::Log::newline(); break;
        case semi_implicit : logging << "bind_stepper(): using semi-implicit stepper" << util::Log::newline(); break;
        case empty :         logging << "bind_stepper(): using empty stepper" << util::Log::newline(); break;
        default :            assert(!"bind_stepper(): unknown stepper type"); break;
        };
    }

    void operator()(auto dt, auto& bodies, const auto& properties, const auto& gravity, const auto& damping,
                    const auto& params, auto& broad_system, auto& narrow_system, auto& contacts, const auto& tag)
    {
        switch (m_type)
        {
        case moreau :
            moreau_time_stepper(dt, bodies, properties, gravity, damping, params, broad_system, narrow_system, contacts,
                                tag);
            break;

        case semi_implicit :
            semi_implicit_time_stepper(dt, bodies, properties, gravity, damping, params, broad_system, narrow_system,
                                       contacts, tag);
            break;

        case empty :
            empty_stepper(dt, bodies, properties, gravity, damping, params, broad_system, narrow_system, contacts, tag);
            break;
        };
    }

private:
    stepper_type m_type;
};

} //namespace prox

// PROX_BIND_STEPPER_H
#endif
