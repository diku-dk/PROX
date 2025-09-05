#ifndef PROX_ENUMS_H
#define PROX_ENUMS_H

namespace prox
{
  /**
   * Different time stepping methods.
   */
enum stepper_type
{
    moreau,
    semi_implicit,
    empty
};

  /**
   * Different sovler types.
   */
enum solver_type
{
    jacobi,
    gauss_seidel
};

  /**
   * Different R-factor strategies for PROX formulation.
   */
enum strategy_type
{
    local_strategy,
    global_strategy,
    blocked_strategy
};

  /**
   * Normal sub solver type.
   * Numerical subroutine used for projecting the normal force onto the normal
   * cone. Different solvers allow for different behaviors of numerical methods.
   */
enum normal_sub_solver_type
{
    nonnegative,
    normal_origin,
    normal_infinity
};

  /**
   * Frictional sub solver type.
   * Numerical subroutine used for projecting the friction force onto the friction
   * cone. Different solvers allow for different shaped friction cones.
   */
enum friction_sub_solver_type
{
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
