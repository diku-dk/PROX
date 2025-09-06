#ifndef PROX_FRICTION_SUB_SOLVER_H
#define PROX_FRICTION_SUB_SOLVER_H

#include "prox_enums.h"
#include "solvers/sub/prox_analytical_ellipsoid.h"
#include "solvers/sub/prox_analytical_sphere.h"
#include "solvers/sub/prox_numerical_ellipsoid.h"
#include "solvers/sub/prox_gjk_ellipsoid.h"
#include "solvers/sub/prox_box_model.h"
#include <stdexcept>

namespace prox
{

template <typename T>
void frictionSolver(friction_sub_solver_type type, const T& z_s, const T& z_t,
                    const T& z_tau, const T& mu_s, const T& mu_t,
                    const T& mu_tau, const T& lambda_n, T& lambda_s,
                    T& lambda_t, T& lambda_tau)
{
    switch (type)
    {
    case analytical_sphere:
        detail::analytical_sphere(z_s, z_t, z_tau, mu_s, lambda_n, lambda_s,
                                  lambda_t, lambda_tau);
        return;

    case analytical_ellipsoid:
        detail::analytical_ellipsoid(z_s, z_t, z_tau, mu_s, mu_tau, lambda_n,
                                     lambda_s, lambda_t, lambda_tau);
        return;

    case numerical_ellipsoid:
        detail::numerical_ellipsoid(z_s, z_t, z_tau, mu_s, mu_t, mu_tau,
                                    lambda_n, lambda_s, lambda_t, lambda_tau);
        return;

    case gjk_ellipsoid:
        detail::gjk_ellipsoid(z_s, z_t, z_tau, mu_s, mu_t, mu_tau, lambda_n,
                              lambda_s, lambda_t, lambda_tau);
        return;

    case box_model:
        detail::box_model(z_s, z_t, z_tau, mu_s, mu_t, mu_tau, lambda_n,
                          lambda_s, lambda_t, lambda_tau);
        return;

    case friction_origin:
        lambda_s = 0;
        lambda_t = 0;
        lambda_tau = 0;
        return;

    case friction_infinity: return;
    default:                throw std::logic_error("Code is broken.");
    }
}

} //namespace prox

// PROX_FRICTION_SUB_SOLVER_H
#endif
