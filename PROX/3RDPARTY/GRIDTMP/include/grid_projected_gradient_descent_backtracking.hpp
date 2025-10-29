#ifndef GRID_PROJECTED_GRADIENT_DESCENT_BACKTRACKING_HPP
#define GRID_PROJECTED_GRADIENT_DESCENT_BACKTRACKING_HPP

#include <grid_grid.h>
#include <grid_local_optimization.hpp>
#include <grid_ccd_frank_wolfe_utils.hpp>
#include <chrono>

namespace grid
{
template <typename T> struct GradientType
{
    T du;
    T dv;
    T dw;
    T dt;
};

template <typename T> void project_to_triangle_barycentric(T& u, T& v, T& w)
{

    //Ensure non-negativity
    u = std::max<T>(0.0, u);
    v = std::max<T>(0.0, v);
    w = std::max<T>(0.0, w);

    //Normalize to sum to 1
    T total = u + v + w;
    if (total > 0)
    {
        u /= total;
        v /= total;
        w /= total;
    }
}

template <typename T>
T project_to_time_interval(const T t, const T tstart, const T tend)
{
    return std::clamp<T>(t, tstart, tend);
}

template <typename T>
EigenVector4<T> compute_gradient(T u, T v, T w, T t,
                                 const RigidBodyInfo<T>& info)
{
    EigenVector3<T> xt = BarycentricInterpolate(u, v, w, t, info);
    //I do not know if we should multiply velocity with t, getVeloctyAtPoint does not use t
    EigenVector3<T> vt = getVelocityAtPoint(info, xt, t) * t;

    //Get SDF gradient at current position
    EigenVector3<T> grad_phi
        = gradientAtProjection(xt, *(info.B_sdf), *(info.B_centerTranslation),
                               *(info.B_centerRotation));

    //Compute vertex positions at time t
    //TODO compute at, bt, ct.
    TriangleAtTimeInfo<T> tri = getTriangleAtTime(t, info);

    //∂φ/∂u
    T du = (grad_phi).dot(tri.A_p0);
    //∂φ/∂v
    T dv = (grad_phi).dot(tri.A_p1);
    //∂φ/∂w
    T dw = (grad_phi).dot(tri.A_p2);
    //∂φ/∂t
    T dt = (grad_phi).dot(vt);

    return EigenVector4<T>(du, dv, dw, dt);
}

template <typename T>
T backtracking_line_search(T tstart, T tend, T u, T v, T w, T t,
                           EigenVector4<T> gradient, EigenVector4<T> direction,
                           const RigidBodyInfo<T>& info, T alpha = 1.0,
                           T beta = 0.5, T c = 1e-4)
{
    EigenVector3<T> xt = BarycentricInterpolate(u, v, w, t, info);
    T current_val
        = valueAtProjection(*(info.B_sdf), xt, *(info.B_centerTranslation),
                            *(info.B_centerRotation));
    //Armijo condition: f(x + αd) ≤ f(x) + cα∇f·d
    T grad_dot_dir = (gradient.dot(direction));

    while (alpha > 1e-10)
    {
        //step
        T u_new = u + alpha * direction.x();
        T v_new = v + alpha * direction.y();
        T w_new = w + alpha * direction.z();
        T t_new = t + alpha * direction.w();

        //Project to triangle
        project_to_triangle_barycentric(u_new, v_new, w_new);
        t_new = project_to_time_interval(t_new, tstart, tend);
        EigenVector3<T> trialPos
            = BarycentricInterpolate(u_new, v_new, w_new, t_new, info);
        T step_val = valueAtProjection(*(info.B_sdf), trialPos,
                                       *(info.B_centerTranslation),
                                       *(info.B_centerRotation));

        if (step_val <= current_val + c * alpha * grad_dot_dir)
        {
            return alpha;
        }

        alpha *= beta;
    }
    //Return smallest step if no better found
    return alpha;
}

template <typename T>
T performProjectedGradientDescent(T tstart, T tend,
                                  const RigidBodyInfo<T>& initialState,
                                  EigenVector3<T>& xtiPoint,
                                  std::vector<T>& minimizerSteps)
{
    minimizerSteps.clear();
    //Initialize with centroid at start time
    T minVel = 0.0000001;
    if ((*(initialState.A_linearVel)).norm() < minVel
        && (*(initialState.A_angularVel)).norm() < minVel)
    {
        return tend;
    }

    EigenVector3<T> p0s = ((initialState.A_p0)).eval();
    EigenVector3<T> p1s = ((initialState.A_p1)).eval();
    EigenVector3<T> p2s = ((initialState.A_p2)).eval();

    EigenVector3<T> vi
        = *(initialState.A_linearVel) - (*(initialState.B_linearVel));
    EigenVector3<T> gradP0 = gradientAtProjection(
        p0s, *(initialState.B_sdf), *(initialState.B_centerTranslation),
        *(initialState.B_centerRotation));
    EigenVector3<T> gradP1 = gradientAtProjection(
        p1s, *(initialState.B_sdf), *(initialState.B_centerTranslation),
        *(initialState.B_centerRotation));
    EigenVector3<T> gradP2 = gradientAtProjection(
        p2s, *(initialState.B_sdf), *(initialState.B_centerTranslation),
        *(initialState.B_centerRotation));
    T p0Min = vi.dot(gradP0);
    T p1Min = vi.dot(gradP1);
    T p2Min = vi.dot(gradP2);
    T u;
    T v;
    T w;
    //I add slight bias such that we will more often select p0 than other vertices
    if (p0Min <= p1Min && p0Min <= p2Min)
    {
        u = 1;
        v = 0;
        w = 0;
    }
    else if (p1Min <= p0Min && p1Min <= p2Min)
    {
        u = 0;
        v = 1;
        w = 0;
    }
    else
    {
        u = 0;
        v = 0;
        w = 1;
    }

    T t = 0.0;
    uint32_t maxIts = 1000;
    T tol = 1e-8;

    for (uint32_t iteration = 0; iteration < maxIts; ++iteration)
    {
        //Compute current gradient
        EigenVector4<T> gradient = compute_gradient(u, v, w, t, initialState);
        /*        EigenVector4<T> gradient
            = {gradientTmp.du, gradientTmp.dv, gradientTmp.dw, gradientTmp.dt};*/

        //Check convergence
        T gradientNorm = (gradient).norm();
        if (gradientNorm < tol) { break; }

        //Descent direction(negative gradient)
        EigenVector4<T> direction = -gradient;

        //Backtracking line search
        auto gss_start = std::chrono::high_resolution_clock::now();
        T alpha = backtracking_line_search(tstart, tend, u, v, w, t, gradient,
                                           direction, initialState);
        auto gss_end = std::chrono::high_resolution_clock::now();
        auto gss_us = std::chrono::duration_cast<std::chrono::nanoseconds>(
                          gss_end - gss_start)
                          .count();
        minimizerSteps.push_back(T(gss_us));

        //Update parameters
        u += alpha * direction.x();
        v += alpha * direction.y();
        w += alpha * direction.z();
        t += alpha * direction.w();

        //Project to constraints
        project_to_triangle_barycentric(u, v, w);
        t = project_to_time_interval(t, tstart, tend);

        //Optional : Check if we've found a collision (φ ≈ 0)
        EigenVector3<T> xti = BarycentricInterpolate(u, v, w, t, initialState);
        EigenVector3<T> vi = getVelocityAtPoint(initialState, xti, t) * t;
        if (std::abs<T>(valueAtProjection(*(initialState.B_sdf), xti,
                                          *(initialState.B_centerTranslation),
                                          *(initialState.B_centerRotation)))
            < tol)
        {
            break;
        }
        xtiPoint = xti;

        bool converged = (iteration < maxIts - 1);
        //return u, v, w, t, converged
        return t;
    }
}
} // namespace grid

#endif // GRID_PROJECTED_GRADIENT_DESCENT_BACKTRACKING_HPP
