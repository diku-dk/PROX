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

#include <solvers/prox_solver.h>
#include <solvers/sub/prox_bind_normal_sub_solver.h>
#include <solvers/sub/prox_bind_friction_sub_solver.h>

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


    //=============EIGEN================0
    Eigen::VectorX<T> qNew; // position vector:               vector7
    Eigen::VectorX<T> uNew; // velocity vector                vector6
    Eigen::VectorX<T> hNew; // external forces:               vector6
    Eigen::VectorX<T> qMNew; // position half step update:    vector7
    Eigen::SparseMatrix<T>
        WNew; // inverse mass matrix:      Diagonal with 6x6 matrices populating only diagonal. Sparse matrix, all other elems 0.
    Eigen::SparseMatrix<T>
        JNew; // Jacobian Matrix:          A matrix where rows = contact, cols = bodies. Sparse matrix. Each contact for a body is a 4x6 matrix. Each non-contact is 0.
    Eigen::SparseMatrix<T>
        WJTNew; // the product of the inverse mass matrix and the transposed Jacobian: Sparse matrix!
    Eigen::VectorX<T> gNew; // correction/stabilization term: vector4
    Eigen::VectorX<T> eNew; // restitution coefficients:      vector4
    Eigen::VectorX<T> wNew; // Current contact velocities:    vector4
    Eigen::VectorX<T>
        WdthNew; // the product of inverse mass matrix, time step and external forces, W*dt*h: vector4
    Eigen::VectorX<T> muNew; // friction coefficients:        vector4
    Eigen::VectorX<T> bNew; // "right hand side vector"       vector4
    Eigen::VectorX<T> lambdaNew; // resulting forces:          vector4
    Eigen::VectorX<T> fcNew; // contact forces:               vector6
    detail::update_body_indices(bodies.begin(), bodies.end());

    get_position_vector_eigen(bodies.begin(), bodies.end(), qNew);

    get_velocity_vector_eigen(bodies.begin(), bodies.end(), uNew);

    if (stepperType == moreau)
    {
        position_update_eigen(qNew, uNew, 0.5f * dt, qMNew);
        //set_position_vector(bodies.begin(), bodies.end(), qM);
        set_position_vector_eigen(bodies.begin(), bodies.end(), qMNew);
    }

    collision_detection(bodies, broad_system, narrow_system, contacts, params);

    unsigned int const number_of_contacts = contacts.size();

    logging << "time_stepper(): Number of contacts = " << number_of_contacts
            << util::Log::newline();

    //EigenSparseVec6<T> hNew;

    get_inverse_mass_matrix_eigen(bodies.begin(), bodies.end(), WNew);

    get_external_forces_vector_eigen(bodies.begin(), bodies.end(), gravity,
                                     damping, hNew);
    //    NCVec6<T> newH = convertEigenToSparse(hNew);

    hNew *= dt;
    WdthNew = WNew * hNew;

    if (number_of_contacts > 0u)
    {

        get_jacobian_matrix_eigen(contacts.begin(), contacts.end(), bodies,
                                  properties, JNew, number_of_contacts);

        if (params.stepper_params().pre_stabilization())
        {

            wNew = JNew * gNew;
            get_pre_stabilization_vector_eigen(contacts.begin(), contacts.end(),
                                               params.stepper_params(), dt,
                                               wNew, gNew, number_of_contacts);
        }
        else
        {
            gNew.resize(number_of_contacts * 4);
            gNew.setZero();
        }

        if (params.stepper_params().bounce_on())
        {
            get_restitution_vector_eigen(contacts.begin(), contacts.end(),
                                         properties, eNew, number_of_contacts);
        }
        else
        {
            //Below *4 because it is vec4s...
            eNew.resize(number_of_contacts * 4);
            eNew.setZero();
        }

        get_friction_coefficient_vector_eigen(contacts.begin(), contacts.end(),
                                              properties, muNew,
                                              number_of_contacts);

        // WJT = M^{-1} J^T
        computeWJT_Eigen(WNew, JNew, WJTNew);

        // b   = (I+E)J u + J W (dt h)
        computeB_Eigen(JNew, WdthNew, uNew, eNew, gNew, bNew);

        run_solver_eigen(JNew, WJTNew, bNew, muNew, lambdaNew,
                         params.solver_params());


        fcNew.resize(WJTNew.rows());
        fcNew.setZero();
        fcNew = WJTNew * lambdaNew;

        // u = u + dt M^{-1} h + fc
        velocity_update_eigen(uNew, WdthNew, fcNew, uNew);
    }
    else
    {
        // u = u + dt M^{-1} h
        velocity_update_eigen(uNew, WdthNew, uNew);
    }

    if (stepperType == moreau)
    {
        position_update_eigen(qMNew, uNew, dt * 0.5f, qNew);
    }
    else if (stepperType == semi_implicit)
    {
        //q^{t+1} = q^t + dt* u^{t+1}
        position_update_eigen(qNew, uNew, dt, qNew);
    }
    //Else empty?

    //set_position_vector(bodies.begin(), bodies.end(), q);
    set_position_vector_eigen(bodies.begin(), bodies.end(), qNew);
    //set_velocity_vector(bodies.begin(), bodies.end(), u);
    set_velocity_vector_eigen(bodies.begin(), bodies.end(), uNew);
    STOP_TIMER("stepper");

    if (params.stepper_params().post_stabilization())
    {
        START_TIMER("stabilization");

        auto newSolverParams = params.solver_params();
        newSolverParams.set_friction_sub_solver(friction_origin);
        newSolverParams.set_normal_sub_solver(nonnegative);

        if (number_of_contacts > 0u)
        {
            get_post_stabilization_vector_eigen(
                contacts.begin(), contacts.end(), params.stepper_params(), gNew,
                number_of_contacts);

            PREFIX("post_");
            run_solver_eigen(JNew, WJTNew, gNew, muNew, lambdaNew,
                             newSolverParams);
            PREFIX("");

            fcNew = WJTNew * lambdaNew;
            position_update_eigen<T>(qNew, fcNew, 1, qNew);

            set_position_vector_eigen(bodies.begin(), bodies.end(), qNew);
        }

        STOP_TIMER("stabilization");
    }
}

template <typename T>
void time_stepper_CCD(T dt, std::vector<RigidBody<T>>& bodies,
                      std::vector<std::vector<Property<T>>> const& properties,
                      Gravity<T> const& gravity, Damping<T> const& damping,
                      Params<T> const& params, broad::System<T>& broad_system,
                      narrow::System<T>& narrow_system,
                      std::vector<ContactPoint<T>>& contacts,
                      std::vector<std::vector<T>>& warmStartBodies)
{
    //auto stepperType = params.stepper_params().stepper();
    stepper_type stepperType = stepper_type::semi_implicit;
    util::Log logging;

    if (stepperType == empty)
    {
        collision_detection(bodies, broad_system, narrow_system, contacts,
                            params);

        logging << "empty_stepper(): Number of contacts = " << contacts.size()
                << util::Log::newline();
        return;
    }

    T tStart = T(0);
    T tEnd = dt;
    T stepDt = dt;

    START_TIMER("Time_Stepper_CCD");
    while (tStart < tEnd)
    {
        START_TIMER("ONE_CCD_ADVANCE");
        collision_detection_only_update_structures(
            bodies, broad_system, narrow_system, contacts, params);

        bool onlyZeroTOI = false;
        START_TIMER("TOI_CCD_FOR_ALL_RIGID_BODIES");
        /*T simulateTo = collision_detection_CCD(bodies, narrow_system, contacts,
                                               tStart, tEnd, onlyZeroTOI);*/
        T simulateTo = collision_detection_CCD_WARM_START(
            bodies, narrow_system, contacts, tStart, tEnd * T(4.0), onlyZeroTOI,
            warmStartBodies);
        simulateTo = std::min<T>(tEnd, simulateTo);
        STOP_TIMER("TOI_CCD_FOR_ALL_RIGID_BODIES");
        std::cerr << "CURRENT simulateTO: " << simulateTo
                  << ", CURRENT tStart, tEnd = (" << tStart << "," << tEnd
                  << ")\n";
        stepDt = simulateTo - tStart;
        if ((std::abs<T>(simulateTo - tEnd) > 0.000001))
        {
            START_TIMER("CCD_DISCRETE_COLLISION_DETECTION_PHASE");
            stepDt = simulateTo - tStart;

            {
                Eigen::VectorX<T> qNew;
                Eigen::VectorX<T> uNew;
                Eigen::VectorX<T> qMNew;
                detail::update_body_indices(bodies.begin(), bodies.end());

                get_position_vector_eigen(bodies.begin(), bodies.end(), qNew);

                get_velocity_vector_eigen(bodies.begin(), bodies.end(), uNew);

                position_update_eigen(qNew, uNew, stepDt, qMNew);
                //set_position_vector(bodies.begin(), bodies.end(), qM);
                set_position_vector_eigen(bodies.begin(), bodies.end(), qMNew);
            }
            tStart = simulateTo;

            collision_detection(bodies, broad_system, narrow_system, contacts,
                                params);
            stepDt = tEnd - tStart;
            //=============EIGEN================0
            Eigen::VectorX<T> qNew; // position vector:               vector7
            Eigen::VectorX<T> uNew; // velocity vector                vector6
            Eigen::VectorX<T> hNew; // external forces:               vector6
            Eigen::VectorX<T> qMNew; // position half step update:    vector7
            Eigen::SparseMatrix<T>
                WNew; // inverse mass matrix:      Diagonal with 6x6 matrices populating only diagonal. Sparse matrix, all other elems 0.
            Eigen::SparseMatrix<T>
                JNew; // Jacobian Matrix:          A matrix where rows = contact, cols = bodies. Sparse matrix. Each contact for a body is a 4x6 matrix. Each non-contact is 0.
            Eigen::SparseMatrix<T>
                WJTNew; // the product of the inverse mass matrix and the transposed Jacobian: Sparse matrix!
            Eigen::VectorX<T> gNew; // correction/stabilization term: vector4
            Eigen::VectorX<T> eNew; // restitution coefficients:      vector4
            Eigen::VectorX<T> wNew; // Current contact velocities:    vector4
            Eigen::VectorX<T>
                WdthNew; // the product of inverse mass matrix, time step and external forces, W*dt*h: vector4
            Eigen::VectorX<T> muNew; // friction coefficients:        vector4
            Eigen::VectorX<T> bNew; // "right hand side vector"       vector4
            Eigen::VectorX<T> lambdaNew; // resulting forces:          vector4
            Eigen::VectorX<T> fcNew; // contact forces:               vector6
            detail::update_body_indices(bodies.begin(), bodies.end());

            get_position_vector_eigen(bodies.begin(), bodies.end(), qNew);

            get_velocity_vector_eigen(bodies.begin(), bodies.end(), uNew);

            unsigned int const number_of_contacts = contacts.size();

            logging << "time_stepper(): Number of contacts = "
                    << number_of_contacts << util::Log::newline();

            //EigenSparseVec6<T> hNew;

            get_inverse_mass_matrix_eigen(bodies.begin(), bodies.end(), WNew);

            get_external_forces_vector_eigen(bodies.begin(), bodies.end(),
                                             gravity, damping, hNew);
            //    NCVec6<T> newH = convertEigenToSparse(hNew);

            hNew *= stepDt;
            WdthNew = WNew * hNew;

            if (number_of_contacts > 0u)
            {
                get_jacobian_matrix_eigen(contacts.begin(), contacts.end(),
                                          bodies, properties, JNew,
                                          number_of_contacts);

                if (params.stepper_params().pre_stabilization())
                {

                    wNew = JNew * gNew;
                    get_pre_stabilization_vector_eigen(
                        contacts.begin(), contacts.end(),
                        params.stepper_params(), stepDt, wNew, gNew,
                        number_of_contacts);
                }
                else
                {
                    gNew.resize(number_of_contacts * 4);
                    gNew.setZero();
                }

                if (params.stepper_params().bounce_on())
                {
                    get_restitution_vector_eigen(contacts.begin(),
                                                 contacts.end(), properties,
                                                 eNew, number_of_contacts);
                }
                else
                {
                    //Below *4 because it is vec4s...
                    eNew.resize(number_of_contacts * 4);
                    eNew.setZero();
                }

                get_friction_coefficient_vector_eigen(
                    contacts.begin(), contacts.end(), properties, muNew,
                    number_of_contacts);

                // WJT = M^{-1} J^T
                computeWJT_Eigen(WNew, JNew, WJTNew);

                // b   = (I+E)J u + J W (dt h)
                computeB_Eigen(JNew, WdthNew, uNew, eNew, gNew, bNew);

                run_solver_eigen(JNew, WJTNew, bNew, muNew, lambdaNew,
                                 params.solver_params());

                fcNew.resize(WJTNew.rows());
                fcNew.setZero();
                fcNew = WJTNew * lambdaNew;

                // u = u + dt M^{-1} h + fc
                velocity_update_eigen(uNew, WdthNew, fcNew, uNew);
            }
            else
            {
                // u = u + dt M^{-1} h
                velocity_update_eigen(uNew, WdthNew, uNew);
            }

            if (stepperType == semi_implicit || true)
            {
                //q^{t+1} = q^t + dt* u^{t+1}
                position_update_eigen(qNew, uNew, stepDt, qNew);
            }
            /*else
            {
                throw std::runtime_error(
                    "EMpty stepper not supported for CCD!");
            }*/
            //Else empty?

            //set_position_vector(bodies.begin(), bodies.end(), q);
            set_position_vector_eigen(bodies.begin(), bodies.end(), qNew);
            //set_velocity_vector(bodies.begin(), bodies.end(), u);
            set_velocity_vector_eigen(bodies.begin(), bodies.end(), uNew);

            if (params.stepper_params().post_stabilization())
            {
                START_TIMER("stabilization");

                auto newSolverParams = params.solver_params();
                newSolverParams.set_friction_sub_solver(friction_origin);
                newSolverParams.set_normal_sub_solver(nonnegative);

                if (number_of_contacts > 0u)
                {
                    get_post_stabilization_vector_eigen(
                        contacts.begin(), contacts.end(),
                        params.stepper_params(), gNew, number_of_contacts);

                    PREFIX("post_");
                    run_solver_eigen(JNew, WJTNew, gNew, muNew, lambdaNew,
                                     newSolverParams);
                    PREFIX("");

                    fcNew = WJTNew * lambdaNew;
                    position_update_eigen<T>(qNew, fcNew, 1, qNew);

                    set_position_vector_eigen(bodies.begin(), bodies.end(),
                                              qNew);
                }

                STOP_TIMER("stabilization");
            }
            tStart += simulateTo;
            STOP_TIMER("CCD_DISCRETE_COLLISION_DETECTION_PHASE");
        }
        else
        {
            if (onlyZeroTOI)
            {
                START_TIMER("CCD_DISCRETE_COLLISION_DETECTION_PHASE_TOI_ZERO");
                //=============EIGEN================0
                Eigen::VectorX<T>
                    qNew; // position vector:               vector7
                Eigen::VectorX<T>
                    uNew; // velocity vector                vector6
                Eigen::VectorX<T>
                    hNew; // external forces:               vector6
                Eigen::VectorX<T>
                    qMNew; // position half step update:    vector7
                Eigen::SparseMatrix<T>
                    WNew; // inverse mass matrix:      Diagonal with 6x6 matrices populating only diagonal. Sparse matrix, all other elems 0.
                Eigen::SparseMatrix<T>
                    JNew; // Jacobian Matrix:          A matrix where rows = contact, cols = bodies. Sparse matrix. Each contact for a body is a 4x6 matrix. Each non-contact is 0.
                Eigen::SparseMatrix<T>
                    WJTNew; // the product of the inverse mass matrix and the transposed Jacobian: Sparse matrix!
                Eigen::VectorX<T>
                    gNew; // correction/stabilization term: vector4
                Eigen::VectorX<T>
                    eNew; // restitution coefficients:      vector4
                Eigen::VectorX<T>
                    wNew; // Current contact velocities:    vector4
                Eigen::VectorX<T>
                    WdthNew; // the product of inverse mass matrix, time step and external forces, W*dt*h: vector4
                Eigen::VectorX<T>
                    muNew; // friction coefficients:        vector4
                Eigen::VectorX<T>
                    bNew; // "right hand side vector"       vector4
                Eigen::VectorX<T>
                    lambdaNew; // resulting forces:          vector4
                Eigen::VectorX<T>
                    fcNew; // contact forces:               vector6
                detail::update_body_indices(bodies.begin(), bodies.end());

                get_position_vector_eigen(bodies.begin(), bodies.end(), qNew);

                get_velocity_vector_eigen(bodies.begin(), bodies.end(), uNew);

                if (stepperType == moreau)
                {
                    position_update_eigen(qNew, uNew, 0.5f * dt, qMNew);
                    //set_position_vector(bodies.begin(), bodies.end(), qM);
                    set_position_vector_eigen(bodies.begin(), bodies.end(),
                                              qMNew);
                }

                collision_detection(bodies, broad_system, narrow_system,
                                    contacts, params);

                unsigned int const number_of_contacts = contacts.size();

                logging << "time_stepper(): Number of contacts = "
                        << number_of_contacts << util::Log::newline();

                //EigenSparseVec6<T> hNew;

                get_inverse_mass_matrix_eigen(bodies.begin(), bodies.end(),
                                              WNew);

                get_external_forces_vector_eigen(bodies.begin(), bodies.end(),
                                                 gravity, damping, hNew);
                //    NCVec6<T> newH = convertEigenToSparse(hNew);

                hNew *= dt;
                WdthNew = WNew * hNew;

                if (number_of_contacts > 0u)
                {

                    get_jacobian_matrix_eigen(contacts.begin(), contacts.end(),
                                              bodies, properties, JNew,
                                              number_of_contacts);

                    if (params.stepper_params().pre_stabilization())
                    {

                        wNew = JNew * gNew;
                        get_pre_stabilization_vector_eigen(
                            contacts.begin(), contacts.end(),
                            params.stepper_params(), dt, wNew, gNew,
                            number_of_contacts);
                    }
                    else
                    {
                        gNew.resize(number_of_contacts * 4);
                        gNew.setZero();
                    }

                    if (params.stepper_params().bounce_on() || true)
                    {
                        get_restitution_vector_eigen(contacts.begin(),
                                                     contacts.end(), properties,
                                                     eNew, number_of_contacts);
                    }
                    else
                    {
                        //Below *4 because it is vec4s...
                        eNew.resize(number_of_contacts * 4);
                        eNew.setZero();
                    }

                    get_friction_coefficient_vector_eigen(
                        contacts.begin(), contacts.end(), properties, muNew,
                        number_of_contacts);

                    // WJT = M^{-1} J^T
                    computeWJT_Eigen(WNew, JNew, WJTNew);

                    // b   = (I+E)J u + J W (dt h)
                    computeB_Eigen(JNew, WdthNew, uNew, eNew, gNew, bNew);

                    run_solver_eigen(JNew, WJTNew, bNew, muNew, lambdaNew,
                                     params.solver_params());

                    fcNew.resize(WJTNew.rows());
                    fcNew.setZero();
                    fcNew = WJTNew * lambdaNew;

                    // u = u + dt M^{-1} h + fc
                    velocity_update_eigen(uNew, WdthNew, fcNew, uNew);
                }
                else
                {
                    // u = u + dt M^{-1} h
                    velocity_update_eigen(uNew, WdthNew, uNew);
                }

                if (stepperType == moreau)
                {
                    position_update_eigen(qMNew, uNew, dt * 0.5f, qNew);
                }
                else if (stepperType == semi_implicit)
                {
                    //q^{t+1} = q^t + dt* u^{t+1}
                    position_update_eigen(qNew, uNew, dt, qNew);
                }
                //Else empty?

                //set_position_vector(bodies.begin(), bodies.end(), q);
                set_position_vector_eigen(bodies.begin(), bodies.end(), qNew);
                //set_velocity_vector(bodies.begin(), bodies.end(), u);
                set_velocity_vector_eigen(bodies.begin(), bodies.end(), uNew);
                STOP_TIMER("stepper");

                if (params.stepper_params().post_stabilization())
                {
                    START_TIMER("stabilization");

                    auto newSolverParams = params.solver_params();
                    newSolverParams.set_friction_sub_solver(friction_origin);
                    newSolverParams.set_normal_sub_solver(nonnegative);

                    if (number_of_contacts > 0u)
                    {
                        get_post_stabilization_vector_eigen(
                            contacts.begin(), contacts.end(),
                            params.stepper_params(), gNew, number_of_contacts);

                        PREFIX("post_");
                        run_solver_eigen(JNew, WJTNew, gNew, muNew, lambdaNew,
                                         newSolverParams);
                        PREFIX("");

                        fcNew = WJTNew * lambdaNew;
                        position_update_eigen<T>(qNew, fcNew, 1, qNew);

                        set_position_vector_eigen(bodies.begin(), bodies.end(),
                                                  qNew);
                    }

                    STOP_TIMER("stabilization");
                }
                STOP_TIMER("CCD_DISCRETE_COLLISION_DETECTION_PHASE_TOI_ZERO");
            }
            else
            {
                START_TIMER("CCD_NO_COLLISION_TOI_IS_END");
                std::cerr << "OTHER CASE!!!!!\n";
                //=============EIGEN================0
                Eigen::VectorX<T>
                    qNew; // position vector:               vector7
                Eigen::VectorX<T>
                    uNew; // velocity vector                vector6
                Eigen::VectorX<T>
                    hNew; // external forces:               vector6
                Eigen::VectorX<T>
                    qMNew; // position half step update:    vector7
                Eigen::SparseMatrix<T>
                    WNew; // inverse mass matrix:      Diagonal with 6x6 matrices populating only diagonal. Sparse matrix, all other elems 0.
                Eigen::SparseMatrix<T>
                    JNew; // Jacobian Matrix:          A matrix where rows = contact, cols = bodies. Sparse matrix. Each contact for a body is a 4x6 matrix. Each non-contact is 0.
                Eigen::SparseMatrix<T>
                    WJTNew; // the product of the inverse mass matrix and the transposed Jacobian: Sparse matrix!
                Eigen::VectorX<T>
                    gNew; // correction/stabilization term: vector4
                Eigen::VectorX<T>
                    eNew; // restitution coefficients:      vector4
                Eigen::VectorX<T>
                    wNew; // Current contact velocities:    vector4
                Eigen::VectorX<T>
                    WdthNew; // the product of inverse mass matrix, time step and external forces, W*dt*h: vector4
                Eigen::VectorX<T>
                    muNew; // friction coefficients:        vector4
                Eigen::VectorX<T>
                    bNew; // "right hand side vector"       vector4
                Eigen::VectorX<T>
                    lambdaNew; // resulting forces:          vector4
                Eigen::VectorX<T>
                    fcNew; // contact forces:               vector6
                detail::update_body_indices(bodies.begin(), bodies.end());

                get_position_vector_eigen(bodies.begin(), bodies.end(), qNew);

                get_velocity_vector_eigen(bodies.begin(), bodies.end(), uNew);

                unsigned int const number_of_contacts = contacts.size();

                logging << "time_stepper(): Number of contacts = "
                        << number_of_contacts << util::Log::newline();

                //EigenSparseVec6<T> hNew;

                get_inverse_mass_matrix_eigen(bodies.begin(), bodies.end(),
                                              WNew);

                get_external_forces_vector_eigen(bodies.begin(), bodies.end(),
                                                 gravity, damping, hNew);
                //    NCVec6<T> newH = convertEigenToSparse(hNew);

                hNew *= stepDt;
                WdthNew = WNew * hNew;

                {
                    // u = u + dt M^{-1} h
                    velocity_update_eigen(uNew, WdthNew, uNew);
                }

                if (stepperType == semi_implicit)
                {
                    //q^{t+1} = q^t + dt* u^{t+1}
                    position_update_eigen(qNew, uNew, stepDt, qNew);
                }
                else
                {
                    throw std::runtime_error(
                        "EMpty stepper not supported for CCD!");
                }
                //Else empty?

                //set_position_vector(bodies.begin(), bodies.end(), q);
                set_position_vector_eigen(bodies.begin(), bodies.end(), qNew);
                //set_velocity_vector(bodies.begin(), bodies.end(), u);
                set_velocity_vector_eigen(bodies.begin(), bodies.end(), uNew);

                //Set start time to end time -- we are finished time stepping!
                tStart = tEnd;
                STOP_TIMER("CCD_NO_COLLISION_TOI_IS_END");
            }
        }
        STOP_TIMER("ONE_CCD_ADVANCE");
    }
    STOP_TIMER("Time_Stepper_CCD");
}

} //namespace prox

// PROX_TIME_STEPPER_H
#endif
