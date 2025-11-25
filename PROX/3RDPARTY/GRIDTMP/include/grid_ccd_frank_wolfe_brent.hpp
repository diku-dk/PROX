#ifndef GRID_CCD_FRANK_WOLFE_BRENT_HPP
#define GRID_CCD_FRANK_WOLFE_BRENT_HPP

#include <grid_grid.h>
#include <grid_local_optimization.hpp>
#include <chrono>

#include <grid_ccd_frank_wolfe_utils.hpp>

namespace grid
{
template <typename T, typename F>
T BrentMinimize_AB(T lstart, T lend, F func,
                   const DistanceAtTimeParams<T>& params,
                   const RigidBodyInfo<T>& info, T tol = 1e-8,
                   int max_iter = 128)
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

    T fx = func(info, params, x);
    T fw = fx;
    T fv = fx;

    T d = (T)0; // current step
    T e = (T)0; // previous step

    for (int iter = 0; iter < max_iter; ++iter)
    {
        T m = (a + c) * (T)0.5; // midpoint of bracket
        T tol1 = tol * std::abs(x) + ZEPS; // absolute tolerance guard
        T tol2 = (T)2 * tol1;

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
            q = (T)2 * (q - r);

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
        T fu = func(info, params, u);

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
T BrentMinimize_ABS(T lstart, T lend, F func,
                    const DistanceAtTimeParams<T>& params,
                    const RigidBodyInfo<T>& info, T tol = 1e-8,
                    int max_iter = 128)
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

    T fx = func(info, params, x);
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
        T fu = func(info, params, u);

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
T BrentMinimize_WORKING(T lstart, T lend, F func,
                        const DistanceAtTimeParams<T>& params,
                        const RigidBodyInfo<T>& info, T tol = 1e-8,
                        int max_iter = 128)
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
        T fa = func(info, params, a);
        T fb = func(info, params, b);
        T fc = func(info, params, c);

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
T BrentMinimize(T lstart, T lend, F func, const DistanceAtTimeParams<T>& params,
                const RigidBodyInfo<T>& info, T tol = 1e-8, int max_iter = 64)
{
    const T golden_ratio = (3 - std::sqrt(T(5))) / T(2);

    T a = lstart;
    T b = lend;
    T x = a + golden_ratio * (b - a);
    T w = x, v = x;

    T fx = func(info, params, x);
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

        T fu = func(info, params, u);

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

    if (func(info, params, x) > 1e-7) return lend;
    return x;
}

template <typename T, typename F>
EigenVector3<T> BrentMinimize(EigenVector3<T> lstart, EigenVector3<T> lend,
                              F func, const DistanceAtPointParams<T>& params,
                              const RigidBodyInfo<T>& info, T tol = T(1e-8),
                              int max_iter = 128)
{
    // Brent's method constants (operating on scalar t in [0,1])
    const T golden_ratio = (3 - std::sqrt(T(5))) / 2; // ~0.38197

    // search on scalar parameter t in [0,1]
    T a = T(0);
    T b = T(1);
    T x = a + golden_ratio * (b - a);
    T w = x, v = x;

    auto point_at
        = [&](T t) -> EigenVector3<T> { return lstart + t * (lend - lstart); };

    T fx = func(info, params, point_at(x));
    T fw = fx, fv = fx;

    T d = T(0), e = T(0);

    for (int iter = 0; iter < max_iter; ++iter)
    {
        T m = (a + b) * T(0.5);
        T tol1 = tol * std::abs(x) + tol;
        T tol2 = T(2) * tol1;

        // Check convergence (on scalar t)
        if (std::abs(x - m) <= tol2 - (b - a) * T(0.5)) { return point_at(x); }

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

            if (q != T(0) && std::abs(p) < std::abs(q * r * T(0.5))
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

        EigenVector3<T> pu = point_at(u);
        T fu = func(info, params, pu);

        // Update brackets (all in scalar t)
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

    // return best found point
    return lstart + x * (lend - lstart);
}

template <typename T>
T FrankWolfeBRENT_BENCHMARK_TIME_STARTING_ITERATE(
    T tstart, T tend, const RigidBodyInfo<T>& initialState,
    DistanceAtTimeParams<T>& distanceAtTimeParams, EigenVector3<T>& xtiPoint,
    std::vector<T>& minimizerSteps)
{
    minimizerSteps.clear();

    //    return FrankWolfeGSSSimple(tstart, tend, initialState);
    T t1 = tstart;
    T ti = t1;
    T end = tend;
    T tip1 = std::numeric_limits<T>::max();
    //For now if there is no velocity, it means our object is stationary. Thus it can never hit the other object. We can thus ignore it.

    float eps = 1e-7;
    size_t maxIterations = 32u;
    //The higher the below iterations, the higher likelihood we will get no false negatives
    size_t hardStopMaxIterations = 1000u;
    size_t its = 0;
    DistanceAtPointParams distanceAtPointParams{.grid = initialState.B_sdf};
    EigenVector3<T> xtip1 = EigenVector3<T>(0, 0, 0);
    EigenVector3<T> xti = EigenVector3<T>(0, 0, 0);

    T minVel = 0.0000001;
    if ((*(initialState.A_linearVel)).norm() < minVel
        && (*(initialState.A_angularVel)).norm() < minVel)
    {
        return tend;
    }
    if (tstart >= tend) { return tend; }

    EigenVector3<T> p0s = ((initialState.A_p0));
    EigenVector3<T> p1s = ((initialState.A_p1));
    EigenVector3<T> p2s = ((initialState.A_p2));

    //For debug purposes to get it out of scope we define it here!
    T phixtip1_2;
    for (size_t i = 0; i < hardStopMaxIterations; ++i)
    {
        its += 1;
        xti = BarycentricInterpolate(distanceAtTimeParams.u,
                                     distanceAtTimeParams.v,
                                     distanceAtTimeParams.w, ti, initialState)
                  .eval();
        //I assume (but only assumption that we call with ti and xti!
        EigenVector3<T> vti = getVelocityAtPoint(initialState, xti, ti);
        T phixti = valueAtProjection(*(initialState.B_sdf), xti,
                                     *(initialState.B_centerTranslation),
                                     *(initialState.B_centerRotation));
        EigenVector3<T> gradPhixti = gradientAtProjection(
            xti, *(initialState.B_sdf), *(initialState.B_centerTranslation),
            *(initialState.B_centerRotation));
        if (phixti <= 0)
        {
            end = std::min<T>(ti, end);
            auto gss_start = std::chrono::high_resolution_clock::now();
            tip1 = BrentMinimize(tstart, ti, UnsignedDistanceAtTime<T>,
                                 distanceAtTimeParams, initialState);
            auto gss_end = std::chrono::high_resolution_clock::now();
            auto gss_us = std::chrono::duration_cast<std::chrono::nanoseconds>(
                              gss_end - gss_start)
                              .count();
            minimizerSteps.push_back(T(gss_us));
        }
        else
        {
            //Compute Direction
            T di = T(-1);
            //Below can be std::copysign not sure. It just says sign in the paper?
            //            if (phixti > 0) { di = -sign(gradPhixti.dot(vi)); }
            if (phixti > 0)
            {
                T val = gradPhixti.dot(vti);
                di = -sign(val);
                //                di = T(1);
            }
            //Direction sign test
            if (di < 0)
            {
                auto gss_start = std::chrono::high_resolution_clock::now();
                tip1 = BrentMinimize(tstart, ti, SignedDistanceAtTime<T>,
                                     distanceAtTimeParams, initialState);
                auto gss_end = std::chrono::high_resolution_clock::now();
                auto gss_us
                    = std::chrono::duration_cast<std::chrono::nanoseconds>(
                          gss_end - gss_start)
                          .count();
                minimizerSteps.push_back(T(gss_us));
            }
            else
            {
                auto gss_start = std::chrono::high_resolution_clock::now();
                tip1 = BrentMinimize(ti, end, SignedDistanceAtTime<T>,
                                     distanceAtTimeParams, initialState);
                auto gss_end = std::chrono::high_resolution_clock::now();
                auto gss_us
                    = std::chrono::duration_cast<std::chrono::nanoseconds>(
                          gss_end - gss_start)
                          .count();
                minimizerSteps.push_back(T(gss_us));
            }
        }
        // Solve spatial sub-problem
        //Note, xtip1 = x_{t_{i+1}}

        xtip1 = BarycentricInterpolate(
                    distanceAtTimeParams.u, distanceAtTimeParams.v,
                    distanceAtTimeParams.w, tip1, initialState)
                    .eval();
        //EigenVector3<T> vtip1 = 0; //TODO compute v_{t_{i+1}}
        EigenVector3<T> vtip1 = getVelocityAtPoint(initialState, xtip1, tip1);
        T phixtip1 = valueAtProjection(*(initialState.B_sdf), xtip1,
                                       *(initialState.B_centerTranslation),
                                       *(initialState.B_centerRotation));
        EigenVector3<T> gradPhixtip1 = gradientAtProjection(
            xtip1, *(initialState.B_sdf), *(initialState.B_centerTranslation),
            *(initialState.B_centerRotation));

        EigenVector3<T> p0_at_ti = getTriangleVertexPosAt(
            *(initialState.A_centerTranslation), *(initialState.A_linearVel),
            *(initialState.A_angularVel), *(initialState.B_centerTranslation),
            *(initialState.B_linearVel), *(initialState.B_angularVel), p0s,
            tip1);
        EigenVector3<T> p1_at_ti = getTriangleVertexPosAt(
            *(initialState.A_centerTranslation), *(initialState.A_linearVel),
            *(initialState.A_angularVel), *(initialState.B_centerTranslation),
            *(initialState.B_linearVel), *(initialState.B_angularVel), p1s,
            tip1);
        EigenVector3<T> p2_at_ti = getTriangleVertexPosAt(
            *(initialState.A_centerTranslation), *(initialState.A_linearVel),
            *(initialState.A_angularVel), *(initialState.B_centerTranslation),
            *(initialState.B_linearVel), *(initialState.B_angularVel), p2s,
            tip1);

        T p0Mins = (p0_at_ti.transpose()).dot(gradPhixtip1);
        T p1Mins = (p1_at_ti.transpose()).dot(gradPhixtip1);
        T p2Mins = (p2_at_ti.transpose()).dot(gradPhixtip1);

        if (phixtip1 <= 0) { end = std::min<T>(tip1, end); }
        EigenVector3<T> si;
        //Computing support vertex pi
        if (p0Mins <= p1Mins && p0Mins <= p2Mins) { si = p0_at_ti; }
        else if (p1Mins <= p2Mins && p1Mins <= p0Mins) { si = p1_at_ti; }
        else { si = p2_at_ti; }
        /*xtip1 = GSSMinimize_WHAT_MODIFIED(T(0), T(1), xtip1, si,
                                          SignedDistanceAtPoint_MODIFIED<T>,
                                          distanceAtPointParams, initialState);*/
        auto gss_start = std::chrono::high_resolution_clock::now();
        xtip1 = GSSMinimize_WHAT(xtip1, si, SignedDistanceAtPoint<T>,
                                 distanceAtPointParams, initialState);
        /*xtip1 = BrentMinimize(xtip1, si, SignedDistanceAtPoint<T>,
                              distanceAtPointParams, initialState);*/
        auto gss_end = std::chrono::high_resolution_clock::now();
        auto gss_us = std::chrono::duration_cast<std::chrono::nanoseconds>(
                          gss_end - gss_start)
                          .count();
        minimizerSteps.push_back(T(gss_us));
        //TODO: Update barycentric coordinates 𝑢, 𝑣, 𝑤 using x®𝑡𝑖+1

        /*computeBarycentricCoordinates(p0_at_ti, p1_at_ti, p2_at_ti, xtip1, u, v,
                                      w);*/
        EigenVector3<T> p0_at_tip1 = getTriangleVertexPosAt(
            *(initialState.A_centerTranslation), *(initialState.A_linearVel),
            *(initialState.A_angularVel), *(initialState.B_centerTranslation),
            *(initialState.B_linearVel), *(initialState.B_angularVel), p0s,
            tip1);
        EigenVector3<T> p1_at_tip1 = getTriangleVertexPosAt(
            *(initialState.A_centerTranslation), *(initialState.A_linearVel),
            *(initialState.A_angularVel), *(initialState.B_centerTranslation),
            *(initialState.B_linearVel), *(initialState.B_angularVel), p1s,
            tip1);
        EigenVector3<T> p2_at_tip1 = getTriangleVertexPosAt(
            *(initialState.A_centerTranslation), *(initialState.A_linearVel),
            *(initialState.A_angularVel), *(initialState.B_centerTranslation),
            *(initialState.B_linearVel), *(initialState.B_angularVel), p2s,
            tip1);
        /*barycentric(p0_at_tip1, p1_at_tip1, p2_at_tip1, xtip1,
                    distanceAtTimeParams.u, distanceAtTimeParams.v,
                    distanceAtTimeParams.w);*/
        computeBarycentricCoordinates(
            p0_at_tip1, p1_at_tip1, p2_at_tip1, xtip1, distanceAtTimeParams.u,
            distanceAtTimeParams.v, distanceAtTimeParams.w);

        projectToTriangle(distanceAtTimeParams.u, distanceAtTimeParams.v,
                          distanceAtTimeParams.w);
        /*distanceAtTimeParams.u = u;
        distanceAtTimeParams.v = v;
        distanceAtTimeParams.w = w;*/

        phixtip1_2 = valueAtProjection(*(initialState.B_sdf), xtip1,
                                       *(initialState.B_centerTranslation),
                                       *(initialState.B_centerRotation));
        if (std::abs(tip1 - ti) <= eps
            && (std::abs(xtip1.x() - xti.x()) <= eps
                && std::abs(xtip1.y() - xti.y()) <= eps
                && std::abs(xtip1.z() - xti.z()) <= eps)
            && phixtip1_2 >= -1e-10 /*&& phixtip1_2 <= eps*/)
        {
            /*if (phixtip1_2 >= eps)
            {
                ti = T(0.01);
                tip1 = T(0.01);
            }*/
            break;
        }
        if (its > maxIterations && phixtip1_2 >= -1e-10) { break; }
        /*if (its > maxIterations)
        {
            std::cerr << "CUrrent it: " << i << " and xtip-value: ("
                      << xtip1.x() << ", " << xtip1.y() << ", " << xtip1.z()
                      << ") and phixtip1_2 =" << phixtip1_2 << "\n";
        }*/
        ti = tip1;
        xti = xtip1;
    }

    xtiPoint = xtip1;
    return std::min<T>(tip1, ti);
}

template <typename T>
T FrankWolfeBRENT_BENCHMARK_TIME_FOUR(T tstart, T tend,
                                      const RigidBodyInfo<T>& initialState,
                                      EigenVector3<T>& xtiPoint,
                                      std::vector<T>& minimizerSteps)
{

    DistanceAtTimeParams distanceAtTimeParams{.u = T(1),
                                              .v = T(0),
                                              .w = T(0),
                                              .grid = (initialState.B_sdf),
                                              .p0 = (initialState.A_p0),
                                              .p1 = (initialState.A_p1),
                                              .p2 = (initialState.A_p2)};
    std::vector<T> tmpMinimizerSteps;

    EigenVector3<T> xtiPointOne;

    T toiOne = FrankWolfeBRENT_BENCHMARK_TIME_STARTING_ITERATE(
        tstart, tend, initialState, distanceAtTimeParams, xtiPointOne,
        tmpMinimizerSteps);
    minimizerSteps.insert(minimizerSteps.end(), tmpMinimizerSteps.begin(),
                          tmpMinimizerSteps.end());
    tmpMinimizerSteps.clear();

    distanceAtTimeParams.u = T(0);
    distanceAtTimeParams.v = T(1);
    distanceAtTimeParams.w = T(0);
    EigenVector3<T> xtiPointTwo;
    T toiTwo = FrankWolfeBRENT_BENCHMARK_TIME_STARTING_ITERATE(
        tstart, tend, initialState, distanceAtTimeParams, xtiPointTwo,
        tmpMinimizerSteps);
    minimizerSteps.insert(minimizerSteps.end(), tmpMinimizerSteps.begin(),
                          tmpMinimizerSteps.end());
    tmpMinimizerSteps.clear();

    distanceAtTimeParams.u = T(0);
    distanceAtTimeParams.v = T(0);
    distanceAtTimeParams.w = T(1);
    EigenVector3<T> xtiPointThree;
    T toiThree = FrankWolfeBRENT_BENCHMARK_TIME_STARTING_ITERATE(
        tstart, tend, initialState, distanceAtTimeParams, xtiPointThree,
        tmpMinimizerSteps);
    minimizerSteps.insert(minimizerSteps.end(), tmpMinimizerSteps.begin(),
                          tmpMinimizerSteps.end());
    tmpMinimizerSteps.clear();

    distanceAtTimeParams.u = T(1.0 / 3.0);
    distanceAtTimeParams.v = T(1.0 / 3.0);
    distanceAtTimeParams.w = T(1.0 / 3.0);
    EigenVector3<T> xtiPointFour;
    T toiFour = FrankWolfeBRENT_BENCHMARK_TIME_STARTING_ITERATE(
        tstart, tend, initialState, distanceAtTimeParams, xtiPointFour,
        tmpMinimizerSteps);
    minimizerSteps.insert(minimizerSteps.end(), tmpMinimizerSteps.begin(),
                          tmpMinimizerSteps.end());
    tmpMinimizerSteps.clear();
    if (toiOne <= toiTwo && toiOne <= toiThree && toiOne <= toiFour
        && valueAtProjection(*(initialState.B_sdf), xtiPointOne,
                             *(initialState.B_centerTranslation),
                             *(initialState.B_centerRotation))
               >= -1e-7)
    {
        xtiPoint = xtiPointOne;
        return toiOne;
    }
    else if (toiTwo <= toiThree && toiTwo <= toiFour
             && valueAtProjection(*(initialState.B_sdf), xtiPointTwo,
                                  *(initialState.B_centerTranslation),
                                  *(initialState.B_centerRotation))
                    >= -1e-7)
    {
        xtiPoint = xtiPointTwo;
        return toiTwo;
    }
    else if (toiThree <= toiFour
             && valueAtProjection(*(initialState.B_sdf), xtiPointThree,
                                  *(initialState.B_centerTranslation),
                                  *(initialState.B_centerRotation))
                    >= -1e-7)
    {
        xtiPoint = xtiPointThree;
        return toiThree;
    }
    else
    {
        xtiPoint = xtiPointFour;
        return toiFour;
    }
}

template <typename T>
T FrankWolfeBRENT_BENCHMARK_TIME_TWO(T tstart, T tend,
                                     const RigidBodyInfo<T>& initialState,
                                     EigenVector3<T>& xtiPoint,
                                     std::vector<T>& minimizerSteps)
{

    EigenVector3<T> p0s = ((initialState.A_p0));
    EigenVector3<T> p1s = ((initialState.A_p1));
    EigenVector3<T> p2s = ((initialState.A_p2));

    //    EigenVector3<T> vi = *(initialState.A_linearVel);
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

    DistanceAtTimeParams distanceAtTimeParams{.u = u,
                                              .v = v,
                                              .w = w,
                                              .grid = (initialState.B_sdf),
                                              .p0 = (initialState.A_p0),
                                              .p1 = (initialState.A_p1),
                                              .p2 = (initialState.A_p2)};
    std::vector<T> tmpMinimizerSteps;

    EigenVector3<T> xtiPointOne;

    T toiOne = FrankWolfeBRENT_BENCHMARK_TIME_STARTING_ITERATE(
        tstart, tend, initialState, distanceAtTimeParams, xtiPointOne,
        tmpMinimizerSteps);
    minimizerSteps.insert(minimizerSteps.end(), tmpMinimizerSteps.begin(),
                          tmpMinimizerSteps.end());
    tmpMinimizerSteps.clear();

    distanceAtTimeParams.u = T(1.0 / 3.0);
    distanceAtTimeParams.v = T(1.0 / 3.0);
    distanceAtTimeParams.w = T(1.0 / 3.0);
    EigenVector3<T> xtiPointFour;
    T toiFour = FrankWolfeBRENT_BENCHMARK_TIME_STARTING_ITERATE(
        tstart, tend, initialState, distanceAtTimeParams, xtiPointFour,
        tmpMinimizerSteps);
    minimizerSteps.insert(minimizerSteps.end(), tmpMinimizerSteps.begin(),
                          tmpMinimizerSteps.end());
    tmpMinimizerSteps.clear();
    if (toiOne <= toiFour
        && valueAtProjection(*(initialState.B_sdf), xtiPointOne,
                             *(initialState.B_centerTranslation),
                             *(initialState.B_centerRotation))
               >= -1e-7)
    {
        xtiPoint = xtiPointOne;
        return toiOne;
    }
    else
    {
        xtiPoint = xtiPointFour;
        return toiFour;
    }
}

template <typename T>
T FrankWolfeBRENT_BENCHMARK_TIME_THREE(T tstart, T tend,
                                       const RigidBodyInfo<T>& initialState,
                                       EigenVector3<T>& xtiPoint,
                                       std::vector<T>& minimizerSteps)
{

    EigenVector3<T> p0s = ((initialState.A_p0));
    EigenVector3<T> p1s = ((initialState.A_p1));
    EigenVector3<T> p2s = ((initialState.A_p2));

    //    EigenVector3<T> vi = *(initialState.A_linearVel);
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
    T u2;
    T v2;
    T w2;
    int c = 0;
    //I add slight bias such that we will more often select p0 than other vertices
    if (p0Min <= p1Min && p0Min <= p2Min)
    {
        u = 1;
        v = 0;
        w = 0;
        if (p1Min <= p2Min)
        {
            u2 = 0;
            v2 = 1;
            w2 = 0;
        }
        else
        {
            u2 = 0;
            v2 = 0;
            w2 = 1;
        }
    }
    else if (p1Min <= p0Min && p1Min <= p2Min)
    {
        u = 0;
        v = 1;
        w = 0;
        if (p0Min <= p2Min)
        {
            u2 = 1;
            v2 = 0;
            w2 = 0;
        }
        else
        {
            u2 = 0;
            v2 = 0;
            w2 = 1;
        }
    }
    else
    {
        if (p0Min <= p1Min)
        {
            u2 = 1;
            v2 = 0;
            w2 = 0;
        }
        else
        {
            u2 = 0;
            v2 = 1;
            w2 = 0;
        }
        u = 0;
        v = 0;
        w = 1;
    }

    DistanceAtTimeParams distanceAtTimeParams{.u = u,
                                              .v = v,
                                              .w = w,
                                              .grid = (initialState.B_sdf),
                                              .p0 = (initialState.A_p0),
                                              .p1 = (initialState.A_p1),
                                              .p2 = (initialState.A_p2)};
    std::vector<T> tmpMinimizerSteps;

    EigenVector3<T> xtiPointOne;

    T toiOne = FrankWolfeBRENT_BENCHMARK_TIME_STARTING_ITERATE(
        tstart, tend, initialState, distanceAtTimeParams, xtiPointOne,
        tmpMinimizerSteps);
    minimizerSteps.insert(minimizerSteps.end(), tmpMinimizerSteps.begin(),
                          tmpMinimizerSteps.end());
    tmpMinimizerSteps.clear();

    EigenVector3<T> xtiPointTwo;
    distanceAtTimeParams.u = u2;
    distanceAtTimeParams.v = v2;
    distanceAtTimeParams.w = w2;
    T toiTwo = FrankWolfeBRENT_BENCHMARK_TIME_STARTING_ITERATE(
        tstart, tend, initialState, distanceAtTimeParams, xtiPointTwo,
        tmpMinimizerSteps);
    minimizerSteps.insert(minimizerSteps.end(), tmpMinimizerSteps.begin(),
                          tmpMinimizerSteps.end());
    tmpMinimizerSteps.clear();

    distanceAtTimeParams.u = T(1.0 / 3.0);
    distanceAtTimeParams.v = T(1.0 / 3.0);
    distanceAtTimeParams.w = T(1.0 / 3.0);
    EigenVector3<T> xtiPointFour;
    T toiFour = FrankWolfeBRENT_BENCHMARK_TIME_STARTING_ITERATE(
        tstart, tend, initialState, distanceAtTimeParams, xtiPointFour,
        tmpMinimizerSteps);
    minimizerSteps.insert(minimizerSteps.end(), tmpMinimizerSteps.begin(),
                          tmpMinimizerSteps.end());
    tmpMinimizerSteps.clear();
    if (toiOne <= toiFour && toiOne <= toiTwo
        && valueAtProjection(*(initialState.B_sdf), xtiPointOne,
                             *(initialState.B_centerTranslation),
                             *(initialState.B_centerRotation))
               >= -1e-7)
    {
        xtiPoint = xtiPointOne;
        return toiOne;
    }
    else if (toiTwo <= toiFour
             && valueAtProjection(*(initialState.B_sdf), xtiPointTwo,
                                  *(initialState.B_centerTranslation),
                                  *(initialState.B_centerRotation))
                    >= -1e-7)
    {
        xtiPoint = xtiPointTwo;
        return toiOne;
    }
    else
    {
        xtiPoint = xtiPointFour;
        return toiFour;
    }
}

template <typename T>
T FrankWolfeBRENT_BENCHMARK_TIME(T tstart, T tend,
                                 const RigidBodyInfo<T>& initialState,
                                 EigenVector3<T>& xtiPoint,
                                 std::vector<T>& minimizerSteps)
{

    EigenVector3<T> p0s = ((initialState.A_p0));
    EigenVector3<T> p1s = ((initialState.A_p1));
    EigenVector3<T> p2s = ((initialState.A_p2));

    //    EigenVector3<T> vi = *(initialState.A_linearVel);
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

    DistanceAtTimeParams distanceAtTimeParams{.u = u,
                                              .v = v,
                                              .w = w,
                                              .grid = (initialState.B_sdf),
                                              .p0 = (initialState.A_p0),
                                              .p1 = (initialState.A_p1),
                                              .p2 = (initialState.A_p2)};
    std::vector<T> tmpMinimizerSteps;

    EigenVector3<T> xtiPointOne;

    T toiOne = FrankWolfeBRENT_BENCHMARK_TIME_STARTING_ITERATE(
        tstart, tend, initialState, distanceAtTimeParams, xtiPointOne,
        tmpMinimizerSteps);
    minimizerSteps.insert(minimizerSteps.end(), tmpMinimizerSteps.begin(),
                          tmpMinimizerSteps.end());
    tmpMinimizerSteps.clear();

    xtiPoint = xtiPointOne;
    return toiOne;
}

template <typename T>
T FrankWolfeBRENT_BENCHMARK_TIME_OLD(T tstart, T tend,
                                     const RigidBodyInfo<T>& initialState,
                                     EigenVector3<T>& xtiPoint,
                                     std::vector<T>& minimizerSteps)
{

    minimizerSteps.clear();

    //    return FrankWolfeGSSSimple(tstart, tend, initialState);
    T t1 = tstart;
    T ti = t1;
    T end = tend;
    T tip1 = std::numeric_limits<T>::max();
    //For now if there is no velocity, it means our object is stationary. Thus it can never hit the other object. We can thus ignore it.
    T minVel = 0.0000001;
    if ((*(initialState.A_linearVel)).norm() < minVel
        && (*(initialState.A_angularVel)).norm() < minVel)
    {
        return tend;
    }
    if (tstart >= tend) { return tend; }

    EigenVector3<T> p0s = ((initialState.A_p0));
    EigenVector3<T> p1s = ((initialState.A_p1));
    EigenVector3<T> p2s = ((initialState.A_p2));

    //    EigenVector3<T> vi = *(initialState.A_linearVel);
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
    /*u = T(1.0 / 3.0);
    v = T(1.0 / 3.0);
    w = T(1.0 / 3.0);*/

    DistanceAtTimeParams distanceAtTimeParams{.u = u,
                                              .v = v,
                                              .w = w,
                                              .grid = (initialState.B_sdf),
                                              .p0 = (initialState.A_p0),
                                              .p1 = (initialState.A_p1),
                                              .p2 = (initialState.A_p2)};

    float eps = 1e-7;
    size_t maxIterations = 32u;
    //The higher the below iterations, the higher likelihood we will get no false negatives
    size_t hardStopMaxIterations = 1000u;
    size_t its = 0;
    DistanceAtPointParams distanceAtPointParams{.grid = initialState.B_sdf};
    EigenVector3<T> xtip1 = EigenVector3<T>(0, 0, 0);
    EigenVector3<T> xti = EigenVector3<T>(0, 0, 0);

    //For debug purposes to get it out of scope we define it here!
    T phixtip1_2;
    for (size_t i = 0; i < hardStopMaxIterations; ++i)
    {
        its += 1;
        xti = BarycentricInterpolate(distanceAtTimeParams.u,
                                     distanceAtTimeParams.v,
                                     distanceAtTimeParams.w, ti, initialState)
                  .eval();
        //I assume (but only assumption that we call with ti and xti!
        EigenVector3<T> vti = getVelocityAtPoint(initialState, xti, ti) * ti;
        T phixti = valueAtProjection(*(initialState.B_sdf), xti,
                                     *(initialState.B_centerTranslation),
                                     *(initialState.B_centerRotation));
        EigenVector3<T> gradPhixti = gradientAtProjection(
            xti, *(initialState.B_sdf), *(initialState.B_centerTranslation),
            *(initialState.B_centerRotation));
        if (phixti <= 0)
        {
            end = std::min<T>(ti, end);
            auto gss_start = std::chrono::high_resolution_clock::now();
            tip1 = BrentMinimize(tstart, ti, UnsignedDistanceAtTime<T>,
                                 distanceAtTimeParams, initialState);
            auto gss_end = std::chrono::high_resolution_clock::now();
            auto gss_us = std::chrono::duration_cast<std::chrono::nanoseconds>(
                              gss_end - gss_start)
                              .count();
            minimizerSteps.push_back(T(gss_us));
        }
        else
        {
            //Compute Direction
            T di = T(-1);
            //Below can be std::copysign not sure. It just says sign in the paper?
            //            if (phixti > 0) { di = -sign(gradPhixti.dot(vi)); }
            if (phixti > 0)
            {
                T val = gradPhixti.dot(vti);
                di = -sign(val);
                //                di = T(1);
            }
            //Direction sign test
            if (di < 0)
            {
                auto gss_start = std::chrono::high_resolution_clock::now();
                tip1 = BrentMinimize(tstart, ti, SignedDistanceAtTime<T>,
                                     distanceAtTimeParams, initialState);
                auto gss_end = std::chrono::high_resolution_clock::now();
                auto gss_us
                    = std::chrono::duration_cast<std::chrono::nanoseconds>(
                          gss_end - gss_start)
                          .count();
                minimizerSteps.push_back(T(gss_us));
            }
            else
            {
                auto gss_start = std::chrono::high_resolution_clock::now();
                tip1 = BrentMinimize(ti, end, SignedDistanceAtTime<T>,
                                     distanceAtTimeParams, initialState);
                auto gss_end = std::chrono::high_resolution_clock::now();
                auto gss_us
                    = std::chrono::duration_cast<std::chrono::nanoseconds>(
                          gss_end - gss_start)
                          .count();
                minimizerSteps.push_back(T(gss_us));
            }
        }
        // Solve spatial sub-problem
        //Note, xtip1 = x_{t_{i+1}}

        xtip1 = BarycentricInterpolate(
                    distanceAtTimeParams.u, distanceAtTimeParams.v,
                    distanceAtTimeParams.w, tip1, initialState)
                    .eval();
        //EigenVector3<T> vtip1 = 0; //TODO compute v_{t_{i+1}}
        EigenVector3<T> vtip1
            = getVelocityAtPoint(initialState, xtip1, tip1) * tip1;
        T phixtip1 = valueAtProjection(*(initialState.B_sdf), xtip1,
                                       *(initialState.B_centerTranslation),
                                       *(initialState.B_centerRotation));
        EigenVector3<T> gradPhixtip1 = gradientAtProjection(
            xtip1, *(initialState.B_sdf), *(initialState.B_centerTranslation),
            *(initialState.B_centerRotation));

        EigenVector3<T> p0_at_ti = getTriangleVertexPosAt(
            *(initialState.A_centerTranslation), *(initialState.A_linearVel),
            *(initialState.A_angularVel), *(initialState.B_centerTranslation),
            *(initialState.B_linearVel), *(initialState.B_angularVel), p0s,
            tip1);
        EigenVector3<T> p1_at_ti = getTriangleVertexPosAt(
            *(initialState.A_centerTranslation), *(initialState.A_linearVel),
            *(initialState.A_angularVel), *(initialState.B_centerTranslation),
            *(initialState.B_linearVel), *(initialState.B_angularVel), p1s,
            tip1);
        EigenVector3<T> p2_at_ti = getTriangleVertexPosAt(
            *(initialState.A_centerTranslation), *(initialState.A_linearVel),
            *(initialState.A_angularVel), *(initialState.B_centerTranslation),
            *(initialState.B_linearVel), *(initialState.B_angularVel), p2s,
            tip1);

        T p0Mins = (p0_at_ti.transpose()).dot(gradPhixtip1);
        T p1Mins = (p1_at_ti.transpose()).dot(gradPhixtip1);
        T p2Mins = (p2_at_ti.transpose()).dot(gradPhixtip1);

        if (phixtip1 <= 0) { end = std::min<T>(tip1, end); }
        EigenVector3<T> si;
        //Computing support vertex pi
        if (p0Mins <= p1Mins && p0Mins <= p2Mins) { si = p0_at_ti; }
        else if (p1Mins <= p2Mins && p1Mins <= p0Mins) { si = p1_at_ti; }
        else { si = p2_at_ti; }
        /*xtip1 = GSSMinimize_WHAT_MODIFIED(T(0), T(1), xtip1, si,
                                          SignedDistanceAtPoint_MODIFIED<T>,
                                          distanceAtPointParams, initialState);*/
        auto gss_start = std::chrono::high_resolution_clock::now();
        /*xtip1 = GSSMinimize_WHAT(xtip1, si, SignedDistanceAtPoint<T>,
                                 distanceAtPointParams, initialState);*/
        xtip1 = BrentMinimize(xtip1, si, SignedDistanceAtPoint<T>,
                              distanceAtPointParams, initialState);
        auto gss_end = std::chrono::high_resolution_clock::now();
        auto gss_us = std::chrono::duration_cast<std::chrono::nanoseconds>(
                          gss_end - gss_start)
                          .count();
        minimizerSteps.push_back(T(gss_us));
        //TODO: Update barycentric coordinates 𝑢, 𝑣, 𝑤 using x®𝑡𝑖+1

        /*computeBarycentricCoordinates(p0_at_ti, p1_at_ti, p2_at_ti, xtip1, u, v,
                                      w);*/
        EigenVector3<T> p0_at_tip1 = getTriangleVertexPosAt(
            *(initialState.A_centerTranslation), *(initialState.A_linearVel),
            *(initialState.A_angularVel), *(initialState.B_centerTranslation),
            *(initialState.B_linearVel), *(initialState.B_angularVel), p0s,
            tip1);
        EigenVector3<T> p1_at_tip1 = getTriangleVertexPosAt(
            *(initialState.A_centerTranslation), *(initialState.A_linearVel),
            *(initialState.A_angularVel), *(initialState.B_centerTranslation),
            *(initialState.B_linearVel), *(initialState.B_angularVel), p1s,
            tip1);
        EigenVector3<T> p2_at_tip1 = getTriangleVertexPosAt(
            *(initialState.A_centerTranslation), *(initialState.A_linearVel),
            *(initialState.A_angularVel), *(initialState.B_centerTranslation),
            *(initialState.B_linearVel), *(initialState.B_angularVel), p2s,
            tip1);
        /*barycentric(p0_at_tip1, p1_at_tip1, p2_at_tip1, xtip1,
                    distanceAtTimeParams.u, distanceAtTimeParams.v,
                    distanceAtTimeParams.w);*/
        computeBarycentricCoordinates(
            p0_at_tip1, p1_at_tip1, p2_at_tip1, xtip1, distanceAtTimeParams.u,
            distanceAtTimeParams.v, distanceAtTimeParams.w);

        projectToTriangle(distanceAtTimeParams.u, distanceAtTimeParams.v,
                          distanceAtTimeParams.w);
        /*distanceAtTimeParams.u = u;
        distanceAtTimeParams.v = v;
        distanceAtTimeParams.w = w;*/

        phixtip1_2 = valueAtProjection(*(initialState.B_sdf), xtip1,
                                       *(initialState.B_centerTranslation),
                                       *(initialState.B_centerRotation));
        if (std::abs(tip1 - ti) <= eps
            && (std::abs(xtip1.x() - xti.x()) <= eps
                && std::abs(xtip1.y() - xti.y()) <= eps
                && std::abs(xtip1.z() - xti.z()) <= eps)
            && phixtip1_2 >= -1e-10 /*&& phixtip1_2 <= eps*/)
        {
            /*if (phixtip1_2 >= eps)
            {
                ti = T(0.01);
                tip1 = T(0.01);
            }*/
            break;
        }
        if (its > maxIterations && phixtip1_2 >= -1e-10) { break; }
        /*if (its > maxIterations)
        {
            std::cerr << "CUrrent it: " << i << " and xtip-value: ("
                      << xtip1.x() << ", " << xtip1.y() << ", " << xtip1.z()
                      << ") and phixtip1_2 =" << phixtip1_2 << "\n";
        }*/
        ti = tip1;
        xti = xtip1;
    }

    xtiPoint = xtip1;
    return tip1;
}
} // namespace grid

#endif
