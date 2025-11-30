#ifndef GRID_SDF_SDF_CCD_GRADIENT_DESCENT_HPP
#define GRID_SDF_SDF_CCD_GRADIENT_DESCENT_HPP

#include "grid_ccd_frank_wolfe_utils.hpp"
#include <grid_grid.h>
#include <grid_sdf_sdf_voxelize.hpp>
#include <chrono>

namespace SDFSDFContact
{
template <typename T> struct currentSDFPose
{
    EigenVector3<T> translation;
    EigenQuaternion<T> rotation;
};

template <typename T> struct SingleRigidBodyInfo
{
    const grid::Grid<T, T>* sdf;
    const EigenVector3<T>* A_linearVel;
    const EigenVector3<T>* A_angularVel;
    const EigenVector3<T>* A_centerTranslation;
    const EigenQuaternion<T>* A_centerRotation;
};

template <typename T>
EigenVector3<T> getVertexPosAt(const EigenVector3<T>& centerTranslation,
                               const EigenVector3<T>& linVel,
                               const EigenVector3<T>& angVel,
                               const EigenVector3<T>& vert, T dt)
{
    //TODO USE EQUATION 14!, DISCARD THIS PELASE!
    //x(t)=R(t)(x−c)+c+v*t
    //And the inverse is x=c+R(t)^{-1}(x(t)-c-v*t)????

    EigenVector3<T> diff = vert - (centerTranslation + linVel * dt);

    return vert + ((linVel + angVel.cross(diff)) * dt);
}

template <typename T>
void getTransformForBody(const EigenVector3<T>& centerTranslation,
                         const EigenVector3<T>& linVel,
                         const EigenVector3<T>& angVel, T dt,
                         EigenVector3<T>& translationVal,
                         EigenQuaternion<T>& quat)
{
    EigenMatrix3<T> rot;
    T eigenNorm = angVel.norm();
    if (std::abs<T>(eigenNorm) > 1e-8)
    {
        const EigenMatrix3<T> skew{
            {static_cast<T>(0),       -angVel.z(),        angVel.y()},
            {       angVel.z(), static_cast<T>(0),       -angVel.x()},
            {      -angVel.y(),        angVel.x(), static_cast<T>(0)}
        };
        rot = EigenMatrix3<T>::Identity()
            + std::sin(eigenNorm * dt) * (skew / eigenNorm)
            + (static_cast<T>(1.0) - std::cos(eigenNorm * dt))
                  * ((skew * skew) / (eigenNorm * eigenNorm));
    }
    else { rot = EigenMatrix3<T>::Identity(); }
    translationVal = linVel * dt;
    EigenQuaternion<T> quaternion(rot);
    quat = quaternion;
}

template <typename T>
EigenVector3<T> getVertexPosAtMat(const EigenVector3<T>& centerTranslation,
                                  const EigenVector3<T>& linVel,
                                  const EigenVector3<T>& angVel,
                                  const EigenVector3<T>& vert, T dt)
{
    EigenMatrix3<T> rot;
    T eigenNorm = angVel.norm();
    if (std::abs<T>(eigenNorm) > 1e-8)
    {
        const EigenMatrix3<T> skew{
            {static_cast<T>(0),       -angVel.z(),        angVel.y()},
            {       angVel.z(), static_cast<T>(0),       -angVel.x()},
            {      -angVel.y(),        angVel.x(), static_cast<T>(0)}
        };
        rot = EigenMatrix3<T>::Identity()
            + std::sin(eigenNorm * dt) * (skew / eigenNorm)
            + (static_cast<T>(1.0) - std::cos(eigenNorm * dt))
                  * ((skew * skew) / (eigenNorm * eigenNorm));
    }
    else { rot = EigenMatrix3<T>::Identity(); }
    return centerTranslation + rot * (vert - centerTranslation) + linVel * dt;
}

template <typename T>
EigenVector3<T> reverseVertexPosAtMat(const EigenVector3<T>& centerTranslation,
                                      const EigenVector3<T>& linVel,
                                      const EigenVector3<T>& angVel,
                                      const EigenVector3<T>& vert, T dt)
{

    EigenMatrix3<T> rot;
    T eigenNorm = angVel.norm();
    //Eigen weirdness???
    /*static_assert(std::is_same_v<float, decltype(std::sin(eigenNorm * dt))>);
    static_assert(std::is_same_v<float, T>);
    static_assert(std::is_same_v<std::complex<float>,
                                 decltype(std::sin<T>(eigenNorm * dt))>);
    static_assert(std::is_same_v<std::complex<float>, T>);
    static_assert(std::is_same_v<float, decltype(std::sin<T>(eigenNorm * dt))>);*/

    if (std::abs<T>(eigenNorm) > 1e-8)
    {
        const EigenMatrix3<T> skew{
            {static_cast<T>(0),       -angVel.z(),        angVel.y()},
            {       angVel.z(), static_cast<T>(0),       -angVel.x()},
            {      -angVel.y(),        angVel.x(), static_cast<T>(0)}
        };
        rot = EigenMatrix3<T>::Identity()
            + std::sin(eigenNorm * dt) * (skew / eigenNorm)
            + (static_cast<T>(1.0) - std::cos(eigenNorm * dt))
                  * ((skew * skew) / (eigenNorm * eigenNorm));
    }
    else { rot = EigenMatrix3<T>::Identity(); }
    return centerTranslation
         + rot.inverse() * (vert - centerTranslation - linVel * dt);
}

template <typename T> T sign(T val)
{
    T res;
    if (val < 0) { res = T(-1); }
    else if (val > 0) { res = T(1); }
    else { res = val; }
    return res;
}

template <typename T, typename F>
T BrentMinimize_EW(T lstart, T lend, F func, const currentSDFPose<T>& poseA,
                   const currentSDFPose<T>& poseB,
                   const SingleRigidBodyInfo<T>& SDFA,
                   const SingleRigidBodyInfo<T>& SDFB,
                   EigenVector3<T> startingPoint, T tol = 1e-10,
                   int max_iter = 64)
{
    // Ensure ordered interval
    T a = lstart;
    T b = lend;
    T c = b;

    const T eps = std::numeric_limits<T>::epsilon();
    // small protective epsilon for tolerance; slightly larger than machine min
    const T ZEPS = std::numeric_limits<T>::min() * (T)1e3;
    T e;
    T d;
    for (int iter = 0; iter < max_iter; ++iter)
    {
        T fa = func(poseA, poseB, SDFA, SDFB, a, startingPoint);
        T fb = func(poseA, poseB, SDFA, SDFB, b, startingPoint);
        T fc = func(poseA, poseB, SDFA, SDFB, c, startingPoint);
        if ((fa > 0 && fb > 0) || fb < 0 && fc < 0)
        {
            c = a;
            d = b - a;
            e = d;
        }
        if (std::abs<T>(fc) < std::abs<T>(fb))
        {
            a = b;
            b = c;
            c = a;
        }
        T xj = (c - b) * T(0.5);
        if (std::abs<T>(xj) < eps || (fb == T(0))) { return b; }
        if ((std::abs<T>(e) >= eps) && (std::abs<T>(fa) > fb))
        {
            T s = fb / fa;
            T p;
            T q;
            if (a == c)
            {
                p = 2 * xj * s;
                q = 1 - s;
            }
            else
            {
                q = fa / fc;
                T r = fb / fc;
                p = s * (2 * xj * q * (q - r) - (b - a) * (r - 1));
                q = (q - 1) * (r - 1) * (s - 1);
            }
            if (p > 0) { q = -q; }
            p = std::abs<T>(p);
            if (2 * p < std::min<T>(3 * xj * q - std::abs<T>(eps * q),
                                    std::abs<T>(e * q)))
            {
                e = d;
                d = p / q;
            }
            else
            {
                d = xj;
                e = d;
            }
        }
        else
        {
            d = xj;
            e = d;
        }
        a = b;
        if (std::abs<T>(d) > eps) { b = b + d; }
        else { b = b + eps * sign(xj); }
    }
    return b;
}

template <typename T, typename F>
T BrentMinimize_BEST(T lstart, T lend, F func, const currentSDFPose<T>& poseA,
                     const currentSDFPose<T>& poseB,
                     const SingleRigidBodyInfo<T>& SDFA,
                     const SingleRigidBodyInfo<T>& SDFB,
                     EigenVector3<T> startingPoint, T tol = 1e-8,
                     int max_iter = 26)
{
    const T golden_ratio
        = (3 - std::sqrt(static_cast<T>(5))) / static_cast<T>(2);

    T a = lstart;
    T b = lend;
    T x = a + golden_ratio * (b - a);
    T w = x, v = x;

    T fx = func(poseA, poseB, SDFA, SDFB, x, startingPoint);
    T fw = fx, fv = fx;

    T d = T(0), e = static_cast<T>(0);

    for (int iter = 0; iter < max_iter; ++iter)
    {
        T m = (a + b) * static_cast<T>(0.5);
        T tol1 = tol * std::abs(x) + tol;
        T tol2 = T(2) * tol1;

        // Check convergence
        if (std::abs(x - m) <= tol2 - (b - a) * static_cast<T>(0.5))
        {
            return x;
        }

        T p = static_cast<T>(0), q = static_cast<T>(0), r = static_cast<T>(0);
        bool do_parabolic = false;

        if (std::abs<T>(e) > tol1)
        {
            // Try parabolic interpolation
            r = (x - w) * (fx - fv);
            q = (x - v) * (fx - fw);
            p = (x - v) * q - (x - w) * r;
            q = T(2) * (q - r);

            if (q > static_cast<T>(0)) { p = -p; }
            else { q = -q; }

            r = e;
            e = d;

            if (std::abs(p) < std::abs(q * r * static_cast<T>(0.5))
                && p > q * (a - x) && p < q * (b - x))
            {
                do_parabolic = true;
                d = p / q;
            }
        }

        if (!do_parabolic)
        {
            // Golden section step
            if (x < m) { e = b - x; }
            else { e = a - x; }
            d = golden_ratio * e;
        }

        T u = x + d;

        // Ensure u is within bounds and not too close to endpoints
        if (u - a < tol2 || b - u < tol2)
        {
            d = (x < m) ? tol1 : -tol1;
            u = x + d;
        }

        T fu = func(poseA, poseB, SDFA, SDFB, u, startingPoint);

        // Update brackets
        if (fu <= fx)
        {
            if (u < x) { b = x; }
            else { a = x; }
            v = w;
            fv = fw;
            w = x;
            fw = fx;
            x = u;
            fx = fu;
        }
        else
        {
            if (u < x) { a = u; }
            else { b = u; }
            if (fu <= fw || w == x)
            {
                v = w;
                fv = fw;
                w = u;
                fw = fu;
            }
            else if (fu <= fv || v == x || v == w)
            {
                v = u;
                fv = fu;
            }
        }
    }

    if (func(poseA, poseB, SDFA, SDFB, x, startingPoint) > 1e-6) return lend;
    return x;
}

template <typename T, typename F>
T BrentMinimize_BAD(T lstart, T lend, F func, const currentSDFPose<T>& poseA,
                    const currentSDFPose<T>& poseB,
                    const SingleRigidBodyInfo<T>& SDFA,
                    const SingleRigidBodyInfo<T>& SDFB,
                    EigenVector3<T> startingPoint, T tol = 1e-8,
                    int max_iter = 1000)
{
    // Ensure ordered interval
    T a = lstart;
    T c = lend;

    const T eps = std::numeric_limits<T>::epsilon();
    // small protective epsilon for tolerance; slightly larger than machine min
    const T ZEPS = std::numeric_limits<T>::min() * (T)1e3;

    // initial guess: midpoint
    T x = (a + c) * (T)0.5;
    T w = x;
    T v = x;

    T fx = func(poseA, poseB, SDFA, SDFB, x, startingPoint);
    T fw = fx;
    T fv = fx;

    T d = (T)0; // current step
    T e = (T)0; // previous step

    for (int iter = 0; iter < max_iter; ++iter)
    {
        T m = (a + c) * (T)0.5; // midpoint of bracket
        T tol1 = tol * std::abs(x) + ZEPS; // absolute tolerance guard
        T tol2 = T(2) * tol1;

        // termination criterion: x is close enough to midpoint
        if (std::abs(x - m) <= tol2 - (c - a) * (T)0.5) { return x; }

        bool used_interpolation = false;
        T p = (T)0, q = (T)0, r = (T)0;
        T u = (T)0; // candidate point to evaluate

        if (std::abs(e) > tol1)
        {
            // Try inverse quadratic interpolation (or secant if appropriate)
            r = (x - w) * (fx - fv);
            q = (x - v) * (fx - fw);
            p = (x - v) * q - (x - w) * r;
            q = T(2) * (q - r);

            if (q != (T)0)
            {
                if (q > (T)0)
                    p = -p;
                else
                    q = -q;

                // Accept interpolation if within bracket and not too large
                if (q != (T)0 && std::abs(p) < std::abs((T)0.5 * q * e)
                    && p > q * (a - x) && p < q * (c - x))
                {
                    d = p / q;
                    u = x + d;

                    // Ensure u isn't too close to boundaries
                    if (u - a < tol2 || c - u < tol2)
                    {
                        // fallback to bisection (but set a small step)
                        d = (m - x);
                        // bisection step: move halfway towards midpoint
                        d = d / (T)2;
                        if (std::abs(d) < tol1) d = (d > (T)0 ? tol1 : -tol1);
                        u = x + d;
                    }

                    used_interpolation = true;
                }
            }
        }

        if (!used_interpolation)
        {
            //BISECTION FALLBACK
            // Move halfway toward the midpoint of [a,c]
            d = (m - x) / (T)2;
            if (std::abs(d) < tol1)
            {
                // ensure we make progress by at least tol1
                d = (m > x) ? tol1 : -tol1;
            }
            u = x + d;
        }

        // Evaluate function at candidate point u
        T fu = func(poseA, poseB, SDFA, SDFB, u, startingPoint);

        // Update bracket and bookkeeping points
        if (fu <= fx)
        {
            // u becomes new best point
            if (u < x)
                c = x;
            else
                a = x;

            // shift v <- w, w <- x, x <- u
            v = w;
            fv = fw;
            w = x;
            fw = fx;
            x = u;
            fx = fu;
        }
        else
        {
            // u is not better than x: shrink bracket toward u
            if (u < x)
                a = u;
            else
                c = u;

            // Update w or v as appropriate
            if (fu <= fw || w == x)
            {
                v = w;
                fv = fw;
                w = u;
                fw = fu;
            }
            else if (fu <= fv || v == x || v == w)
            {
                v = u;
                fv = fu;
            }
        }

        // Prepare for next iteration: keep last step magnitude in e
        e = d;
    }

    // max_iter reached: return best found x
    return x;
}

template <typename T, typename F>
T BrentMinimize(T lstart, T lend, F func, const currentSDFPose<T>& poseA,
                const currentSDFPose<T>& poseB,
                const SingleRigidBodyInfo<T>& SDFA,
                const SingleRigidBodyInfo<T>& SDFB,
                EigenVector3<T> startingPoint, int max_iter = 128)
{
    const T tol = std::sqrt(std::numeric_limits<T>::epsilon());

    // Quick initial sampling (fewer points for speed)
    const int quick_samples = 5;
    T best_x = (lstart + lend) / 2.0;
    T best_fx = func(poseA, poseB, SDFA, SDFB, best_x, startingPoint);

    for (int i = 0; i <= quick_samples; ++i)
    {
        T x = lstart + (lend - lstart) * i / quick_samples;
        T fx = func(poseA, poseB, SDFA, SDFB, x, startingPoint);
        if (fx < best_fx)
        {
            best_x = x;
            best_fx = fx;
        }
    }

    // Use the best point and its neighbors for parabolic interpolation
    T a = std::max(lstart, best_x - (lend - lstart) / quick_samples);
    T b = best_x;
    T c = std::min(lend, best_x + (lend - lstart) / quick_samples);

    T fa = func(poseA, poseB, SDFA, SDFB, a, startingPoint);
    T fb = best_fx;
    T fc = func(poseA, poseB, SDFA, SDFB, c, startingPoint);

    // Main parabolic interpolation loop
    for (int iter = 0; iter < max_iter; ++iter)
    {
        // Try parabolic interpolation
        T denom = (a - b) * (a - c) * (b - c);
        if (std::abs(denom) < tol)
        {
            // Points are too close, use bisection
            T x_new = (a + c) / 2.0;
        }
        else
        {
            T A = (c * (fb - fa) + b * (fa - fc) + a * (fc - fb)) / denom;
            T B = (c * c * (fa - fb) + b * b * (fc - fa) + a * a * (fb - fc))
                / denom;

            if (A > tol)
            { // Convex parabola - has a minimum
                T x_new = -B / (2.0 * A);

                // Check if the new point is reasonable
                if (x_new > a && x_new < c && std::abs(x_new - b) > tol)
                {
                    T f_new
                        = func(poseA, poseB, SDFA, SDFB, x_new, startingPoint);

                    // Update the triplet
                    if (f_new < fb)
                    {
                        // New point is better, shift points
                        if (x_new > b)
                        {
                            a = b;
                            fa = fb;
                            b = x_new;
                            fb = f_new;
                        }
                        else
                        {
                            c = b;
                            fc = fb;
                            b = x_new;
                            fb = f_new;
                        }
                    }
                    else
                    {
                        // New point is worse, contract the interval
                        if (x_new < b)
                        {
                            a = x_new;
                            fa = f_new;
                        }
                        else
                        {
                            c = x_new;
                            fc = f_new;
                        }
                    }
                    continue;
                }
            }
        }

        // If parabolic interpolation failed, use a fallback strategy
        T step = (c - a) / 4.0;
        T candidates[3] = {a + step, b, c - step};
        T best_candidate = b;
        T best_fval = fb;

        for (int i = 0; i < 3; ++i)
        {
            if (i == 1) continue; // Skip b since we already know its value

            T fval
                = func(poseA, poseB, SDFA, SDFB, candidates[i], startingPoint);
            if (fval < best_fval)
            {
                best_candidate = candidates[i];
                best_fval = fval;
            }
        }

        // Update the interval around the best candidate
        T new_step = (c - a) / 6.0;
        a = std::max(lstart, best_candidate - new_step);
        c = std::min(lend, best_candidate + new_step);
        b = best_candidate;
        fa = func(poseA, poseB, SDFA, SDFB, a, startingPoint);
        fb = best_fval;
        fc = func(poseA, poseB, SDFA, SDFB, c, startingPoint);

        // Check convergence
        if (c - a < tol * (std::abs(b) + 1)) { break; }
    }

    return b;
}

template <typename T, typename F>
T BrentMinimize_SLOW_BUT_ACCURATE(T lstart, T lend, F func,
                                  const currentSDFPose<T>& poseA,
                                  const currentSDFPose<T>& poseB,
                                  const SingleRigidBodyInfo<T>& SDFA,
                                  const SingleRigidBodyInfo<T>& SDFB,
                                  EigenVector3<T> startingPoint,
                                  int max_iter = 1000)
{

    // Use a hybrid approach: initial global sampling + local pattern search
    const int num_initial_samples = 35;
    const int pattern_search_iterations = 60;
    const T convergence_tol = std::sqrt(std::numeric_limits<T>::epsilon());

    // Phase 1: Initial global sampling to find promising regions
    std::vector<T> sample_points;
    std::vector<T> sample_values;

    T best_x = lstart;
    T best_fx = func(poseA, poseB, SDFA, SDFB, lstart, startingPoint);

    // Sample the function at multiple points
    for (int i = 0; i <= num_initial_samples; ++i)
    {
        T x = lstart + (lend - lstart) * i / num_initial_samples;
        T fx = func(poseA, poseB, SDFA, SDFB, x, startingPoint);
        sample_points.push_back(x);
        sample_values.push_back(fx);

        if (fx < best_fx)
        {
            best_x = x;
            best_fx = fx;
        }
    }

    // Phase 2: Pattern search from the best points found
    auto pattern_search = [&](T start_x)
    {
        T x = start_x;
        T fx = func(poseA, poseB, SDFA, SDFB, start_x, startingPoint);
        T step_size = (lend - lstart) / 8.0;

        for (int iter = 0;
             iter < pattern_search_iterations && step_size > convergence_tol;
             ++iter)
        {
            bool improved = false;

            // Explore positive direction
            T x_plus = x + step_size;
            if (x_plus <= lend)
            {
                T f_plus
                    = func(poseA, poseB, SDFA, SDFB, x_plus, startingPoint);
                if (f_plus < fx)
                {
                    x = x_plus;
                    fx = f_plus;
                    improved = true;
                }
            }

            // Explore negative direction
            T x_minus = x - step_size;
            if (x_minus >= lstart)
            {
                T f_minus
                    = func(poseA, poseB, SDFA, SDFB, x_minus, startingPoint);
                if (f_minus < fx)
                {
                    x = x_minus;
                    fx = f_minus;
                    improved = true;
                }
            }

            // If no improvement, reduce step size
            if (!improved) { step_size *= 0.5; }
        }

        return std::make_pair(x, fx);
    };

    // Run pattern search from the top 3 initial samples
    std::vector<std::pair<T, T>> candidates;
    for (int i = 0; i < std::min(3, num_initial_samples + 1); ++i)
    {
        // Find the next best starting point that we haven't used yet
        T start_x = sample_points[0];
        T start_fx = sample_values[0];
        int best_idx = 0;

        for (size_t j = 1; j < sample_points.size(); ++j)
        {
            bool already_used = false;
            for (const auto& cand : candidates)
            {
                if (std::abs(sample_points[j] - cand.first) < convergence_tol)
                {
                    already_used = true;
                    break;
                }
            }

            if (!already_used && sample_values[j] < start_fx)
            {
                start_x = sample_points[j];
                start_fx = sample_values[j];
                best_idx = j;
            }
        }

        auto result = pattern_search(start_x);
        candidates.push_back(result);

        // Update global best
        if (result.second < best_fx)
        {
            best_x = result.first;
            best_fx = result.second;
        }
    }

    // Phase 3: Final refinement using quadratic interpolation around the best point
    T final_x = best_x;
    T final_fx = best_fx;

    for (int refine_iter = 0; refine_iter < 5; ++refine_iter)
    {
        T h = (lend - lstart) * std::pow(0.5, refine_iter + 2);

        T x1 = std::max(lstart, final_x - h);
        T x2 = final_x;
        T x3 = std::min(lend, final_x + h);

        T f1 = func(poseA, poseB, SDFA, SDFB, x1, startingPoint);
        T f2 = final_fx;
        T f3 = func(poseA, poseB, SDFA, SDFB, x3, startingPoint);

        // Try quadratic interpolation if we have three distinct points
        if (std::abs(x1 - x2) > convergence_tol
            && std::abs(x2 - x3) > convergence_tol
            && std::abs(x1 - x3) > convergence_tol)
        {

            T denom = (x1 - x2) * (x1 - x3) * (x2 - x3);
            if (std::abs(denom) > convergence_tol)
            {
                T a = (x3 * (f2 - f1) + x2 * (f1 - f3) + x1 * (f3 - f2))
                    / denom;
                T b = (x3 * x3 * (f1 - f2) + x2 * x2 * (f3 - f1)
                       + x1 * x1 * (f2 - f3))
                    / denom;

                if (a > 0)
                { // Convex parabola
                    T x_min = -b / (2 * a);
                    if (x_min >= x1 && x_min <= x3)
                    {
                        T f_min = func(poseA, poseB, SDFA, SDFB, x_min,
                                       startingPoint);
                        if (f_min < final_fx)
                        {
                            final_x = x_min;
                            final_fx = f_min;
                            continue;
                        }
                    }
                }
            }
        }

        // If quadratic interpolation didn't help, do a local pattern search
        T step = h / 4.0;
        bool improved = false;

        for (int dir = -1; dir <= 1; dir += 2)
        {
            T test_x = final_x + dir * step;
            if (test_x >= lstart && test_x <= lend)
            {
                T test_fx
                    = func(poseA, poseB, SDFA, SDFB, test_x, startingPoint);
                if (test_fx < final_fx)
                {
                    final_x = test_x;
                    final_fx = test_fx;
                    improved = true;
                }
            }
        }

        if (!improved) { break; }
    }

    //    if (func(poseA, poseB, SDFA, SDFB, final_x, startingPoint) > 1e-6) return lend;
    return final_x;
}

template <typename T, typename F>
T BrentMinimize_OK(T lstart, T lend, F func, const currentSDFPose<T>& poseA,
                   const currentSDFPose<T>& poseB,
                   const SingleRigidBodyInfo<T>& SDFA,
                   const SingleRigidBodyInfo<T>& SDFB,
                   EigenVector3<T> startingPoint, T tol = 1e-8,
                   int max_iter = 1000)
{
    const T golden_ratio = (3.50 - std::sqrt(T(6))) / T(2);

    T a = lstart;
    T b = lend;
    T x = golden_ratio * (b - a);
    T w = x, v = x;

    T fx = func(poseA, poseB, SDFA, SDFB, x, startingPoint);
    T fw = fx, fv = fx;

    T d = T(0), e = T(0);

    for (int iter = 0; iter < max_iter; ++iter)
    {
        T m = (a + b) * T(0.5);
        T tol1 = tol * std::abs(x) + tol;
        T tol2 = T(2) * tol1;

        // Check convergence
        if (std::abs(x - m) <= tol2 - (b - a) * T(0.5)) { return x; }

        T p = T(0), q = T(0), r = T(0);
        bool do_parabolic = false;

        if (std::abs(e) > tol1)
        {
            // Try parabolic interpolation
            r = (x - w) * (fx - fv);
            q = (x - v) * (fx - fw);
            p = (x - v) * q - (x - w) * r;
            q = T(2) * (q - r);

            if (q > T(0)) { p = -p; }
            else { q = -q; }

            r = e;
            e = d;

            if (std::abs(p) < std::abs(q * r * T(0.5)) && p > q * (a - x)
                && p < q * (b - x))
            {
                do_parabolic = true;
                d = p / q;
            }
        }

        if (!do_parabolic)
        {
            //BISECTION FALLBACK
            // Move halfway toward the midpoint of [a,c]
            d = (m - x) / (T)2;
            if (std::abs(d) < tol1)
            {
                // ensure we make progress by at least tol1
                d = (m > x) ? tol1 : -tol1;
            }
        }

        T u = x + d;

        // Ensure u is within bounds and not too close to endpoints
        if (u - a < tol2 || b - u < tol2)
        {
            d = (x < m) ? tol1 : -tol1;
            u = x + d;
        }

        T fu = func(poseA, poseB, SDFA, SDFB, u, startingPoint);

        // Update brackets
        if (fu <= fx)
        {
            if (u < x) { b = x; }
            else { a = x; }
            v = w;
            fv = fw;
            w = x;
            fw = fx;
            x = u;
            fx = fu;
        }
        else
        {
            if (u < x) { a = u; }
            else { b = u; }
            if (fu <= fw || w == x)
            {
                v = w;
                fv = fw;
                w = u;
                fw = fu;
            }
            else if (fu <= fv || v == x || v == w)
            {
                v = u;
                fv = fu;
            }
        }
    }

    if (func(poseA, poseB, SDFA, SDFB, x, startingPoint) > 1e-7) return lend;
    return x;
}

template <typename T>
EigenVector3<T> projectToSDFSurfaceLocal(EigenVector3<T> pos,
                                         const SingleRigidBodyInfo<T>& sdfInfo)
{
    return pos
         - grid::valueAtProjection(*(sdfInfo.sdf), pos,
                                   *(sdfInfo.A_centerTranslation),
                                   *(sdfInfo.A_centerRotation))
               * grid::gradientAtProjection(pos, *(sdfInfo.sdf),
                                            *(sdfInfo.A_centerTranslation),
                                            *(sdfInfo.A_centerRotation));
}

//Apply the standard SDF transformation and then the time transformation
template <typename T>
EigenVector3<T> gradientAtProjectionForB(const EigenVector3<T>& samplePoint,
                                         const grid::Grid<T, T>& sdf,
                                         const EigenVector3<T>& translationSDF,
                                         const EigenQuaternion<T>& rotationSDF,
                                         const currentSDFPose<T> pose)
{
    EigenVector3<T> localSamplePoint
        = pose.rotation.inverse()
        * ((rotationSDF.inverse() * (samplePoint - translationSDF))
           - pose.translation);
    /*    contactPoint
        = transformRotation * contactPoint + transformTranslation;*/
    return (pose.rotation * rotationSDF)
         * (grid::computeGradient_Working(localSamplePoint, sdf));
}

//Apply the standard SDF transformation and then the time transformation
template <typename T>
T valueAtProjectionForB(const EigenVector3<T>& samplePoint,
                        const grid::Grid<T, T>& sdf,
                        const EigenVector3<T>& translationSDF,
                        const EigenQuaternion<T>& rotationSDF,
                        const currentSDFPose<T> pose)
{
    EigenVector3<T> localSamplePoint
        = pose.rotation.inverse()
        * ((rotationSDF.inverse() * (samplePoint - translationSDF))
           - pose.translation);
    return (grid::value_at_2<T, T>(sdf, localSamplePoint));
}

template <typename T>
T SignedDistance(const currentSDFPose<T>& poseA_dead,
                 const currentSDFPose<T>& poseB_dead,
                 const SingleRigidBodyInfo<T>& SDFA,
                 const SingleRigidBodyInfo<T>& SDFB, T ti,
                 EigenVector3<T> point)
{
    //Not sure if this is the actual correct way to do it!
    //We always have point sticks to SDF A, thus these should always move towards eachother
    //Thus we don't have to simulate the WHOLE SDF -- only the point with the SDF A's linear and angular velocity
    EigenVector3<T> newPos
        = getVertexPosAtMat(*(SDFA.A_centerTranslation), *(SDFA.A_linearVel),
                            *(SDFA.A_angularVel), point, ti);
    //But of course newPos is in time ti, thus this is completely in world coordinates. Henceforth,
    // we need to "simulate" sdf pose at time t, by getting its translation and rotation at time t.
    // then we can use this pose to transform from world -> local:
    EigenVector3<T> translationB;
    EigenQuaternion<T> rotationB;
    getTransformForBody(*(SDFB.A_centerTranslation), *(SDFB.A_linearVel),
                        *(SDFB.A_angularVel), ti, translationB, rotationB);
    currentSDFPose<T> poseB = {translationB, rotationB};

    /*    return std::abs<T>(grid::valueAtProjection(
        SDFB, newPos, *(SDFB.A_centerTranslation), *(SDFB.A_centerRotation)));*/
    return std::abs<T>(valueAtProjectionForB(newPos, *(SDFB.sdf),
                                             *(SDFB.A_centerTranslation),
                                             *(SDFB.A_centerRotation), poseB));
}

template <typename T>
T getSDFSDFTOISingleVoxelCCD(const EigenVector3<T>& position,
                             const SingleRigidBodyInfo<T>& SDFA,
                             const SingleRigidBodyInfo<T>& SDFB, T tstart,
                             T tend, EigenVector3<T>& outContactPoint)
{
    //TODO: If velocity is zero for this SDF combination, we can return because nothing can move in time. Maybe return -1 for nothing? Otherwise return 1 for end

    //First our pos is in local coordinates in the sdf with absolutely no transformations
    // not even those applied to the SDF. Thus we get the pos in world coordinates at time t=0:
    EigenVector3<T> pos
        = *(SDFA.A_centerRotation) * position + *(SDFA.A_centerTranslation);
    T ti = tstart;
    T tip1 = -1.0;
    EigenVector3<T> xip1;
    EigenVector3<T> x_ti;
    uint32_t maxIterations = 10000;
    T stepSizeAlpha = 0.0001;
    T eps = 1e-8;
    for (uint32_t i = 0; i < maxIterations; ++i)
    {
        //Assume pos is in world coordinates but not transformed. So it is in world coordinates compared to
        // SDF_A. So basically pos is always the point in world coordinates at time t=0.

        //Let x always be in local coordiantes. Then we need no transformation for x
        // in SDF A. x should only be transformed to B when we do gradient computation with B
        EigenVector3<T> translationB;
        EigenQuaternion<T> rotationB;
        EigenVector3<T> translationA;
        EigenQuaternion<T> rotationA;
        //Transform pos from time t=0 to ti.
        //x_ti is in WORLD coordinates!
        x_ti = getVertexPosAtMat(*(SDFA.A_centerTranslation),
                                 *(SDFA.A_linearVel), *(SDFA.A_angularVel), pos,
                                 ti);
        //We do the same as before, but ONLY get the rotation and translation needed
        // to transform SDF B into time ti, such that we can make accurate quries for
        // x_ti.
        getTransformForBody(*(SDFB.A_centerTranslation), *(SDFB.A_linearVel),
                            *(SDFB.A_angularVel), ti, translationB, rotationB);
        getTransformForBody(*(SDFA.A_centerTranslation), *(SDFA.A_linearVel),
                            *(SDFA.A_angularVel), ti, translationA, rotationA);
        currentSDFPose<T> poseB = {translationB, rotationB};
        currentSDFPose<T> poseA = {translationA, rotationA};

        //WRONG BELOW?
        //Uhhh we never update rotation. Let us first check translation works, then test rotation.
        //Yeah change below to projection based on new translation!

        if (valueAtProjectionForB(x_ti, *(SDFB.sdf),
                                  *(SDFB.A_centerTranslation),
                                  *(SDFB.A_centerRotation), poseB)
            < static_cast<T>(0))

        {
            tend = std::min<T>(tend, ti);
            tip1 = BrentMinimize<T>(tstart, ti, SignedDistance<T>, poseA, poseB,
                                    SDFA, SDFB, pos);
        }
        else
        {
            //            auto start = std::chrono::high_resolution_clock::now();
            tip1 = BrentMinimize_SLOW_BUT_ACCURATE<T>(ti, tend * 0.5,
                                                      SignedDistance<T>, poseA,
                                                      poseB, SDFA, SDFB, pos);
            T tip2 = BrentMinimize_SLOW_BUT_ACCURATE<T>(
                tend * 0.5, tend, SignedDistance<T>, poseA, poseB, SDFA, SDFB,
                pos);
            if (std::abs<T>(tip1 - tend * 0.5) < eps) { tip1 = tip2; }
            else { tip1 = std::min<T>(tip1, tip2); }
            return tip1;

            /*            auto end = std::chrono::high_resolution_clock::now();
            auto duration
                = std::chrono::duration_cast<std::chrono::nanoseconds>(end
                                                                       - start);
            uint64_t ns_duration = duration.count();
            std::cerr << ns_duration << ", ";
            return tip1;*/
        }

        //The normal is quite simple to compute. We know pos is in world coordinates, but at time t=0.
        //So we only need to transform it with the SDF's local coordinates!
        //NO ABOVE IS INCORRECT: IT SHOULD BE AT TIME t=ti!!!
        /*EigenVector3<T> normA
            = grid::gradientAtProjection(pos, *(SDFA.sdf),
                                         *(SDFA.A_centerTranslation),
                                         *(SDFA.A_centerRotation))
                  .normalized();*/
        EigenVector3<T> normA
            = gradientAtProjectionForB(x_ti, *(SDFA.sdf),
                                       *(SDFA.A_centerTranslation),
                                       *(SDFA.A_centerRotation), poseA)
                  .normalized();

        //Now for GRADB it is different! First we assume x_ti is pos at time ti.
        //That means real world coordinates, so now we should apply the transformations in
        //reverse that made pos move to time t_i. So we computed the rotation and translation
        // of SDFB in time t_i and move the x_ti -> local SDF
        EigenVector3<T> gradB = gradientAtProjectionForB(
            x_ti, *(SDFB.sdf), *(SDFB.A_centerTranslation),
            *(SDFB.A_centerRotation), poseB);
        EigenVector3<T> gradientDir = (gradB - (gradB.dot(normA)) * normA);
        xip1 = x_ti - stepSizeAlpha * gradientDir;
        //Check for convergence...
        //Maybe check tip?
        if (std::abs<T>(xip1.norm() - x_ti.norm()) < eps
            || std::abs<T>(tip1 - ti) <= eps)
        {
            break;
        }

        //Now traverse back to SDF start pose, such that xtip now lies in the
        // SDFs pose at t=0!
        xip1 = reverseVertexPosAtMat(*(SDFA.A_centerTranslation),
                                     *(SDFA.A_linearVel), *(SDFA.A_angularVel),
                                     xip1, ti);
        //After transform, remember to project back to local coordinates!
        xip1 = projectToSDFSurfaceLocal(xip1, SDFA);
        //Now set our new search start point to pos!
        pos = xip1;
    }
    outContactPoint = xip1;
    return tip1;
}

template <typename T>
T getSDFSDFTOISingleVoxel(const SDFVoxel<T>& SDFBVoxel,
                          const SingleRigidBodyInfo<T>& SDFA,
                          const SingleRigidBodyInfo<T>& SDFB, T tstart, T tend,
                          std::vector<EigenVector3<T>>& outContactPoints)
{
    T minTOI = std::numeric_limits<T>::max();
    for (size_t i = 0; i < SDFBVoxel.selected.size(); ++i)
    {
        EigenVector3<T> selectedPoint = SDFBVoxel.selected[i].pos;
        EigenVector3<T> outContactPoint;
        T toi = getSDFSDFTOISingleVoxelCCD(selectedPoint, SDFA, SDFB, tstart,
                                           tend, outContactPoint);
        if (std::abs<T>(toi - minTOI) <= 1e-9)
        {
            outContactPoints.push_back(outContactPoint);
        }
        else if (toi <= minTOI)
        {
            outContactPoints.clear();
            outContactPoints.push_back(outContactPoint);
            minTOI = toi;
        }
    }
    return minTOI;
}

template <typename T>
T getSDFSDFTOI(const std::vector<SDFVoxel<T>>& SDFAVoxel,
               const std::vector<SDFVoxel<T>>& SDFBVoxel,
               const SingleRigidBodyInfo<T>& SDFA,
               const SingleRigidBodyInfo<T>& SDFB, T tstart, T tend,
               std::vector<EigenVector3<T>>& outContactPoints)
{

    T minTOI = std::numeric_limits<T>::max();
    //Compute all voxels of SDF A
    for (size_t i = 0; i < SDFAVoxel.size(); ++i)
    {
        std::vector<EigenVector3<T>> newOutContactPoints;
        //Importantly to get correct projection first loop is SDFA, SDFB, second should swap!
        T toi = getSDFSDFTOISingleVoxel(SDFAVoxel[i], SDFA, SDFB, tstart, tend,
                                        newOutContactPoints);
        if (minTOI >= toi)
        {
            minTOI = toi;
            outContactPoints = std::move(newOutContactPoints);
        };
    }

    //Compute all voxels of SDF B
    for (size_t i = 0; i < SDFBVoxel.size(); ++i)
    {
        std::vector<EigenVector3<T>> newOutContactPoints;
        //Importantly to get correct projection first loop is SDFA, SDFB, second should swap!
        T toi = getSDFSDFTOISingleVoxel(SDFBVoxel[i], SDFB, SDFA, tstart, tend,
                                        newOutContactPoints);
        if (minTOI >= toi)
        {
            minTOI = toi;
            outContactPoints = std::move(newOutContactPoints);
        }
    }
    return minTOI;
}
} // namespace SDFSDFContact

#endif // GRID_SDF_SDF_CCD_GRADIENT_DESCENT_HPP
