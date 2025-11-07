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
                           T rho = 0.5, T c = 1e-4)
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
        T u_new = u;
        T v_new = v;
        T w_new = w;
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

        alpha *= rho;
    }
    //Return smallest step if no better found
    return alpha;
}

template <typename T>
T performProjectedGradientDescentOld(T tstart, T tend,
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
        //        if (gradientNorm < tol) { break; }

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

        //Check if we've found a collision (φ ≈ 0)
        EigenVector3<T> xti = BarycentricInterpolate(u, v, w, t, initialState);
        EigenVector3<T> vi = getVelocityAtPoint(initialState, xti, t) * t;
        T distTo = std::abs<T>(valueAtProjection(
            *(initialState.B_sdf), xti, *(initialState.B_centerTranslation),
            *(initialState.B_centerRotation)));
        if (distTo < tol) { break; }
        xtiPoint = xti;

        bool converged = (iteration < maxIts - 1);
        //return u, v, w, t, converged
    }
    return t;
}

template <typename T>
T performProjectedGradientDescent_TMP(T tstart, T tend,
                                      const RigidBodyInfo<T>& initialState,
                                      EigenVector3<T>& contactPoint,
                                      std::vector<T>& minimizerSteps)
{
    T minVel = 0.0000001;
    if ((*(initialState.A_linearVel)).norm() < minVel
        && (*(initialState.A_angularVel)).norm() < minVel)
    {
        return tend;
    }

    T u, v, w;

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

    T t = tstart; // Start searching from beginning of interval
    uint32_t maxIts = 100;
    T tol = 1e-6;

    // Store the best candidate found so far
    T best_t = tend;
    T best_sdf = std::numeric_limits<T>::max();
    EigenVector3<T> best_point;
    bool found_penetration = false;

    for (uint32_t iteration = 0; iteration < maxIts; ++iteration)
    {
        EigenVector3<T> xt = BarycentricInterpolate(u, v, w, t, initialState);
        T sdf_val = valueAtProjection(*(initialState.B_sdf), xt,
                                      *(initialState.B_centerTranslation),
                                      *(initialState.B_centerRotation));

        // CRITICAL: Check if we've found a collision (φ ≈ 0 from positive side)
        if (sdf_val >= -tol && sdf_val <= tol && !found_penetration)
        {
            contactPoint = xt;
            return t; // Found first time of impact!
        }

        // If we're penetrating, we've gone too far
        if (sdf_val < -tol)
        {
            found_penetration = true;
            // Move backward in time to find the actual TOI
            best_t = t; // This is our new upper bound
            t = (tstart + t) * 0.5; // Binary search backward
            continue;
        }

        // Store the best candidate so far (closest to zero from positive side)
        if (sdf_val > 0 && sdf_val < best_sdf)
        {
            best_sdf = sdf_val;
            best_t = t;
            best_point = xt;
        }

        // Compute gradient and take a step
        EigenVector4<T> gradient = compute_gradient(u, v, w, t, initialState);

        // For CCD, we want to move toward φ=0, not necessarily minimize φ
        // So we use the negative gradient but with careful step sizing
        EigenVector4<T> direction = -gradient;
        auto gss_start = std::chrono::high_resolution_clock::now();
        T alpha
            = backtracking_line_search(tstart, tend, u, v, w, t, gradient,
                                       direction, initialState, 0.1, 0.5, 1e-4);
        auto gss_end = std::chrono::high_resolution_clock::now();
        auto gss_us = std::chrono::duration_cast<std::chrono::nanoseconds>(
                          gss_end - gss_start)
                          .count();
        minimizerSteps.push_back(T(gss_us));
        // Update with smaller, more conservative steps for TOI finding
        u += 0.5 * alpha * direction.x();
        v += 0.5 * alpha * direction.y();
        w += 0.5 * alpha * direction.z();
        t += 0.5 * alpha * direction.w();

        // Project to constraints
        project_to_triangle_barycentric(u, v, w);
        //t = project_to_time_interval(t, tstart, tend);

        // Check convergence
        //if (gradient.norm() < tol) { break; }
    }

    // If we found a good candidate close to the surface, use it
    if (best_sdf < 0.01)
    { // Within 1cm of surface
        contactPoint = best_point;
        return best_t;
    }

    return best_t; // No collision found
}

template <typename T> void ProjectToSimplex(T& u, T& v, T& w)
{
    // Put coordinates into a vector and sort descending
    std::vector<T> vals = {u, v, w};
    std::vector<T> mu = vals;
    std::sort(mu.begin(), mu.end(), std::greater<double>());

    // Find rho = largest index j such that mu[j] > (sum_{k<=j} mu[k] - 1) / (j+1)
    T cumulative = 0.0;
    int rho = -1;
    T cumulative_at_rho = 0.0;
    for (int j = 0; j < 3; ++j)
    {
        cumulative += mu[j];
        T t = (cumulative - 1.0)
            / (j + 1); // candidate theta for first j+1 entries
        if (mu[j] > t)
        {
            rho = j;
            cumulative_at_rho = cumulative;
        }
    }

    // Degenerate case: no positive rho found -> assign centroid (1/3,1/3,1/3)
    if (rho == -1)
    {
        u = T(1.0 / 3.0);
        v = T(1.0 / 3.0);
        w = T(1.0 / 3.0);
        return;
    }

    // Compute theta and project
    T theta = (cumulative_at_rho - 1.0) / (rho + 1);
    Eigen::Vector3d proj;
    proj[0] = std::max(u - theta, 0.0);
    proj[1] = std::max(v - theta, 0.0);
    proj[2] = std::max(w - theta, 0.0);

    // Normalize to guard against tiny FP error so the sum is exactly 1
    double s = proj.sum();
    if (s > 0.0)
        proj /= s;
    else
        proj = {T(1.0 / 3.0), T(1.0 / 3.0), T(1.0 / 3.0)};

    u = proj[0];
    v = proj[1];
    w = proj[2];
    return;
}

template <typename T>
T performProjectedGradientDescent(T tstart, T tend,
                                  const RigidBodyInfo<T>& initialState,
                                  EigenVector3<T>& contactPoint,
                                  std::vector<T>& minimizerSteps)
{
    T minVel = 0.0000001;
    if ((*(initialState.A_linearVel)).norm() < minVel
        && (*(initialState.A_angularVel)).norm() < minVel)
    {
        return tend;
    }

    T u, v, w;

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

    T t = tstart; // Start searching from beginning of interval
    uint32_t maxIts = 1000;
    T tol = 1e-6;

    // Store the best candidate found so far
    T best_t = tend;
    T best_sdf = std::numeric_limits<T>::max();
    EigenVector3<T> best_point;
    bool found_penetration = false;

    EigenVector4<T> initialSolution = {u, v, w, t};
    EigenVector4<T> initialSolutionLast
        = {std::numeric_limits<T>::max(), std::numeric_limits<T>::max(),
           std::numeric_limits<T>::max(), std::numeric_limits<T>::max()};
    T phival;
    for (uint32_t iteration = 0; iteration < maxIts; ++iteration)
    {
        EigenVector3<T> xt = BarycentricInterpolate(
            initialSolution.x(), initialSolution.y(), initialSolution.z(),
            initialSolution.w(), initialState);
        phival = valueAtProjection(*(initialState.B_sdf), xt,
                                   *(initialState.B_centerTranslation),
                                   *(initialState.B_centerRotation));

        if (phival <= 1e-6 && phival > -1e-6) { break; }
        EigenVector3<T> vt = getVelocityAtPoint(initialState, xt, t);

        TriangleAtTimeInfo<T> tri = getTriangleAtTime(t, initialState);

        EigenVector3<T> gradPhi = gradientAtProjection(
            xt, *(initialState.B_sdf), *(initialState.B_centerTranslation),
            *(initialState.B_centerRotation));
        T du = dot(gradPhi, tri.A_p0);
        T dv = dot(gradPhi, tri.A_p1);
        T dw = dot(gradPhi, tri.A_p2);
        T dt = dot(gradPhi, vt);

        EigenVector4<T> gradient = {du, dv, dw, dt};
        //Whenever we get an initial TOI, we know our gradient does not need to search further, as this TOI below is the upper bound!
        /*if (phival < 0)
        {
            tend = std::min<T>(tend, initialSolution.w());
            dt = -std::abs<T>(dt);
            gradient = {du, dv, dw, dt};
        }*/

        //Backtracking LINE search!
        T alpha = 1.0;
        //T alpha = T(1.0) / (std::abs<T>(dt) / tend);
        T rho = 0.5; // Reduction factor for backtracking (0 < rho < 1)
        T c_armijo
            = 1e-4; // Armijo condition constant (sufficient decrease parameter)
        T tol = 1e-6; // Convergence tolerance
        EigenVector4<T> c = initialSolution;
        EigenVector4<T> p = -gradient;
        T deltaFk = (gradient.transpose()).dot(p);
        bool sufficientDecrease = true;
        T alpha_min = 1e-10;
        alpha = alpha_min;
        /*        while (alpha > 1e-10)
        {
            EigenVector4<T> cnew = initialSolution + alpha * p;
            ProjectToSimplex(cnew.x(), cnew.y(), cnew.z());
            cnew.w() = std::clamp<T>(cnew.w(), tstart, tend);
            EigenVector3<T> xnew = BarycentricInterpolate(
                cnew.x(), cnew.y(), cnew.z(), cnew.w(), initialState);
            T phiNew = valueAtProjection(*(initialState.B_sdf), xnew,
                                         *(initialState.B_centerTranslation),
                                         *(initialState.B_centerRotation));
            if (phiNew <= phival + c_armijo * alpha * deltaFk) { break; }
            else { alpha *= rho; }
        }*/
        bool found_any = false;
        int max_iters = 20;
        T alpha_max = 1.0;
        T prev_alpha = 0; // last accepted alpha
        T eta = 1.1;
        /*if (deltaFk >= 0)
        {
            // Not a descent direction: do not search
            alpha = 0;
        }
        else*/
        {
            auto gss_start = std::chrono::high_resolution_clock::now();
            // First make sure alpha_min actually changes the projected point:
            for (int iter = 0; iter < max_iters && alpha <= alpha_max; ++iter)
            {
                EigenVector4<T> cnew = initialSolution + alpha * p;
                ProjectToSimplex(cnew.x(), cnew.y(), cnew.z());
                cnew.w() = std::clamp<T>(cnew.w(), tstart, tend);

                EigenVector3<T> xnew = BarycentricInterpolate(
                    cnew.x(), cnew.y(), cnew.z(), cnew.w(), initialState);
                T phiNew
                    = valueAtProjection(*(initialState.B_sdf), xnew,
                                        *(initialState.B_centerTranslation),
                                        *(initialState.B_centerRotation));

                if (phiNew <= phival + c_armijo * alpha * deltaFk)
                {
                    // condition holds: record and try a larger alpha
                    prev_alpha = alpha;
                    found_any = true;
                    alpha *= eta;
                }
                else
                {
                    // condition failed: we have bracket (prev_alpha, alpha)
                    alpha = prev_alpha;
                    break;
                }
            }

            if (found_any && false)
            {
                // Optionally refine using bisection between prev_alpha and min(alpha,alpha_max)
                T left = prev_alpha;
                T right = std::min(alpha, alpha_max);
                for (int k = 0; k < 10; ++k)
                { // a few bisection steps
                    T mid = (left + right) * 0.5;
                    EigenVector4<T> cnew = initialSolution + mid * p;
                    ProjectToSimplex(cnew.x(), cnew.y(), cnew.z());
                    cnew.w() = std::clamp<T>(cnew.w(), tstart, tend);
                    if ((cnew - initialSolution).norm() == 0)
                    {
                        left = mid;
                        continue;
                    }
                    EigenVector3<T> xnew = BarycentricInterpolate(
                        cnew.x(), cnew.y(), cnew.z(), cnew.w(), initialState);
                    T phiNew
                        = valueAtProjection(*(initialState.B_sdf), xnew,
                                            *(initialState.B_centerTranslation),
                                            *(initialState.B_centerRotation));
                    if (phiNew <= phival + c_armijo * mid * deltaFk)
                    {
                        left = mid;
                    }
                    else { right = mid; }
                }
                alpha = left; // accept refined alpha
            }
            auto gss_end = std::chrono::high_resolution_clock::now();
            auto gss_us = std::chrono::duration_cast<std::chrono::nanoseconds>(
                              gss_end - gss_start)
                              .count();
            minimizerSteps.push_back(T(gss_us));
        }
        //alpha = (T(1.0) / (std::abs<T>(dt) / tend)) * 0.0001;
        //alpha = 1e-6;

        EigenVector4<T> tmpSolution = initialSolution - alpha * gradient;
        ProjectToSimplex(tmpSolution.x(), tmpSolution.y(), tmpSolution.z());
        tmpSolution.w() = std::clamp<T>(tmpSolution.w(), tstart, tend);
        initialSolution = tmpSolution;

        //if ((initialSolution - initialSolutionLast).norm() < 1e-7) { break; }
        initialSolutionLast = initialSolution;
    }
    if (std::abs<T>(phival) < tol) { ; }
    contactPoint = BarycentricInterpolate(
        initialSolution.x(), initialSolution.y(), initialSolution.z(),
        initialSolution.w(), initialState);
    return initialSolution.w();

    //    return tend; // No collision found
}

template <typename T>
T performProjectedGradientDescent_AHA2(T tstart, T tend,
                                       const RigidBodyInfo<T>& initialState,
                                       EigenVector3<T>& contactPoint,
                                       std::vector<T>& minimizerSteps)
{
    const int num_starts = 5; // Try 5 starting points
    T best_t = tend;
    EigenVector3<T> best_contact;
    T best_sdf = std::numeric_limits<T>::max();

    for (int i = 0; i < num_starts; ++i)
    {
        T endt = ((num_starts + 1) - (i + 1)) * (tend / T(num_starts));
        EigenVector3<T> current_contact;
        T current_t = performProjectedGradientDescent_AHA(
            tstart, endt, initialState, current_contact, minimizerSteps);

        // Evaluate the solution
        EigenVector3<T> xt = current_contact;
        T current_sdf = valueAtProjection(*(initialState.B_sdf), xt,
                                          *(initialState.B_centerTranslation),
                                          *(initialState.B_centerRotation));

        std::cerr << "CURRDISTTOSDF: " << current_sdf << " and " << current_t
                  << "\n";
        if (((std::abs(current_sdf) < std::abs(best_sdf)))
            && current_t < best_t)
        {
            best_t = current_t;
            best_contact = current_contact;
            best_sdf = current_sdf;
        }
    }

    contactPoint = best_contact;
    return best_t;
}
} // namespace grid

#endif // GRID_PROJECTED_GRADIENT_DESCENT_BACKTRACKING_HPP
