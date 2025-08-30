#ifndef PROX_ENUMS_H
#define PROX_ENUMS_H

namespace prox
{
  /**
   * Different time stepping methods.
   */
using stepper_type = enum { moreau, semi_implicit, empty };

  /**
   * Different sovler types.
   */
using solver_type = enum { jacobi, gauss_seidel };

  /**
   * Different R-factor strategies for PROX formulation.
   */
using strategy_type = enum { local_strategy, global_strategy, blocked_strategy };

  /**
   * Normal sub solver type.
   * Numerical subroutine used for projecting the normal force onto the normal
   * cone. Different solvers allow for different behaviors of numerical methods.
   */
using normal_sub_solver_type = enum { nonnegative, normal_origin, normal_infinity };

  /**
   * Frictional sub solver type.
   * Numerical subroutine used for projecting the friction force onto the friction
   * cone. Different solvers allow for different shaped friction cones.
   */
using friction_sub_solver_type = enum {
    analytical_sphere,
    analytical_ellipsoid,
    numerical_ellipsoid,
    gjk_ellipsoid,
    box_model,
    friction_origin,
    friction_infinity
};

} // namespace prox

// PROX_ENUMS_H
#endif
