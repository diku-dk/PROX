#ifndef PROX_TIME_STEPPER_H
#define PROX_TIME_STEPPER_H

#include <prox_rigid_body.h>
#include <prox_contact_point.h>

#include <prox_get_mass_matrix.h>
#include <prox_get_inverse_mass_matrix.h>
#include <prox_get_jacobian_matrix.h>
#include <prox_get_external_forces_vector.h>
#include <prox_get_position_vector.h>
#include <prox_get_velocity_vector.h>
#include <prox_get_pre_stabilization_vector.h>
#include <prox_get_post_stabilization_vector.h>
#include <prox_get_restitution_vector.h>
#include <prox_get_friction_coefficient_vector.h>

#include <prox_set_position_vector.h>
#include <prox_set_velocity_vector.h>

#include <prox_position_update.h>
#include <prox_velocity_update.h>
#include <prox_collision_detection.h>
#include <prox_update_body_indices.h>

#include <prox_params.h>
#include <prox_math_policy.h>

#include <solvers/prox_bind_solver.h>
#include <solvers/sub/prox_bind_normal_sub_solver.h>
#include <solvers/sub/prox_bind_friction_sub_solver.h>
#include <solvers/strategies/prox_bind_R_strategy.h>

#include <util_profiling.h>
#include <broad.h>
#include <narrow.h>

#include <vector>

namespace prox
{

template <typename T>
void time_stepper(T dt, std::vector<RigidBody<T>>& bodies,
                  std::vector<std::vector<Property<T>>> const& properties,
                  Gravity<T> const& gravity, Damping<T> const& damping,
                  Params<T> const& params, broad::System<T>& broad_system,
                  narrow::System<T>& narrow_system,
                  std::vector<ContactPoint<T>>& contacts)
{
    auto stepperType = params.stepper_params().stepper();
    util::Log logging;

    START_TIMER("stepper");

    if (stepperType == empty)
    {
        collision_detection(bodies, broad_system, narrow_system, contacts,
                            params);

        logging << "empty_stepper(): Number of contacts = " << contacts.size()
                << util::Log::newline();

        STOP_TIMER("stepper");
        return;
    }

    NCVec7<T> q;      // position vector
    NCVec7<T> qM;     // position half step update
    NCVec6<T> u;      // velocity vector
    DiagonalMatrix<6, T> W;      // inverse mass matrix.
        // 2009-08-13 Kenny code review: Optimization replace diagonal6x6_type
        // with diagonal_mass_type, maybe wait to optimize until all it working
    NCVec6<T> h; // external forces
    NCVec6<T>
        Wdth; // the product of inverse mass matrix, time step and external forces, W*dt*h
    NCVec6<T> fc; // contact forces
    NCVec4<T> g; // correction/stabilization term
    NCVec4<T> e; // restitution coefficients
    NCVec4<T> mu; // friction coefficients
    NCVec4<T> lambda; // resulting forces
    NCVec4<T> b; // "right hand side vector"
    CRMatrix<4, 6, T> J; // Jacobian matrix
    CRMatrix<6, 4, T>
        WJT; // the product of the inverse mass matrix and the transposed Jacobian
    NCVec4<T> w; // Current contact velocities

    detail::update_body_indices(bodies.begin(), bodies.end());

    get_position_vector(bodies.begin(), bodies.end(), q);
    get_velocity_vector(bodies.begin(), bodies.end(), u);

    if (stepperType == moreau)
    {
        position_update(q, u, .5f * dt, qM);
        set_position_vector(bodies.begin(), bodies.end(), qM);
    }

    collision_detection(bodies, broad_system, narrow_system, contacts, params);

    unsigned int const number_of_contacts = contacts.size();

    logging << "thyme_stepper(): Number of contacts = " << number_of_contacts
            << util::Log::newline();

    get_inverse_mass_matrix(bodies.begin(), bodies.end(), W);
    get_external_forces_vector(bodies.begin(), bodies.end(), gravity, damping,
                               h);

    sparse::prod(dt, h);
    sparse::prod(W, h, Wdth); // Wdth = dt M^{-1} f_ext

    if (number_of_contacts > 0u)
    {
        get_jacobian_matrix(contacts.begin(), contacts.end(), bodies,
                            properties, J, number_of_contacts);

        if (params.stepper_params().pre_stabilization())
        {
            sparse::prod(J, u, w);

            get_pre_stabilization_vector(contacts.begin(), contacts.end(),
                                         params.stepper_params(), dt, w, g,
                                         number_of_contacts);
        }
        else
        {
            g.resize(number_of_contacts);
            g.clear_data();
        }

        if (params.stepper_params().bounce_on())
        {
            get_restitution_vector(contacts.begin(), contacts.end(), properties,
                                   e, number_of_contacts);
        }
        else
        {
            e.resize(number_of_contacts);
            e.clear_data();
        }

        get_friction_coefficient_vector(contacts.begin(), contacts.end(),
                                        properties, mu, number_of_contacts);

        computeWJT(W, J, WJT); // WJT = M^{-1} J^T

        computeB(J, Wdth, u, e, g, b); // b   = (I+E)J u + J W (dt h)

        run_solver(J, WJT, b, mu, lambda, params.solver_params());

        fc.resize(WJT.nrows());

        sparse::prod(WJT, lambda, fc, true);     // fc = M^{-1}*J^T*lambda

        velocity_update(u, Wdth, fc, u); // u = u + dt M^{-1} h + fc
    }
    else
    {
        velocity_update(u, Wdth, u); // u = u + dt M^{-1} h
    }

    if (stepperType == moreau) { position_update(qM, u, dt * .5f, q); }
    else if (stepperType == semi_implicit) { position_update(q, u, dt, q); }

    set_position_vector(bodies.begin(), bodies.end(), q);
    set_velocity_vector(bodies.begin(), bodies.end(), u);

    STOP_TIMER("stepper");

    if (params.stepper_params().post_stabilization())
    {
        START_TIMER("stabilization");

        auto newSolverParams = params.solver_params();
        newSolverParams.set_friction_sub_solver(friction_origin);
        newSolverParams.set_normal_sub_solver(nonnegative);

        if (number_of_contacts > 0u)
        {
            get_post_stabilization_vector(contacts.begin(), contacts.end(),
                                          params.stepper_params(), g,
                                          number_of_contacts);

            PREFIX("post_");
            run_solver(J, WJT, g, mu, lambda, newSolverParams);
            PREFIX("");

            sparse::prod(WJT, lambda, fc, true);

            position_update<T>(q, fc, 1, q);

            set_position_vector(bodies.begin(), bodies.end(), q);
        }

        STOP_TIMER("stabilization");
    }
}

} //namespace prox

// PROX_TIME_STEPPER_H
#endif
