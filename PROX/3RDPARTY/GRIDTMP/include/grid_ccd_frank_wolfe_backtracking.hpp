#ifndef GRID_CCD_FRANK_WOLFE_BACKTRACKNG_HPP
#define GRID_CCD_FRANK_WOLFE_BACKTRACKNG_HPP

#include <grid_grid.h>
#include <grid_local_optimization.hpp>
#include <chrono>

#include <grid_ccd_frank_wolfe_utils.hpp>

namespace grid
{
template <typename T, typename F>
T armijoTemporalSearch(F f, const DistanceAtTimeParams<T>& params,
                       const RigidBodyInfo<T>& info, T tLow, T tHigh,
                       T tCurrent, T fCurrent, T dfCurrent)
{
    T alphaInit = 1.0;
    T c = 0.0009;
    T rho = 0.5;
    int maxArmijoIters = 40;

    //Determine search direction based on original algorithm logic
    T direction = (dfCurrent < 0) ? 1.0 : -1.0;

    T alpha = alphaInit;
    T bestT = tCurrent;

    for (int i = 0; i < maxArmijoIters; ++i)
    {
        //We use relative scaling, because we have an upper and lower bound we must respect!
        T tCandidate = tCurrent + direction * alpha * (tHigh - tLow);

        //Bound check!
        if (tCandidate < tLow) { tCandidate = tLow; }
        if (tCandidate > tHigh) { tCandidate = tHigh; }

        T fCandidate = f(info, params, tCandidate);

        //Armijo condition for minimization: f(t + αd) ≤ f(t) + c * α * ∇f·d
        //Note: ∇f·d = dfCurrent * direction
        T armijoBound = fCurrent + c * alpha * dfCurrent * direction;

        if (fCandidate <= armijoBound)
        {
            bestT = tCandidate;
            break;
        }

        alpha *= rho;

        //Early termination
        if (alpha < 1e-10) break;
    }

    return bestT;
}

template <typename T, typename F>
EigenVector3<T> armijoSpatialSearch(F f, const DistanceAtPointParams<T>& params,
                                    const RigidBodyInfo<T>& info,
                                    const EigenVector3<T>& x_current,
                                    const EigenVector3<T>& support_vertex)
{
    T alphaInit = 1.0;
    T c = 0.0009;
    T rho = 0.5;
    int maxArmijoIters = 40;
    EigenVector3<T> direction = support_vertex - x_current;
    T alpha = alphaInit;

    T f0 = valueAtProjection(*(params.grid), x_current,
                             *(info.B_centerTranslation),
                             *(info.B_centerRotation));
    EigenVector3<T> grad = gradientAtProjection(x_current, *(params.grid),
                                                *(info.B_centerTranslation),
                                                *(info.B_centerRotation));
    T dfDirection = grad.dot(direction);

    for (int i = 0; i < maxArmijoIters; ++i)
    {
        Eigen::Vector3<T> xCandidate = x_current + alpha * direction;

        // Project back to triangle if needed (assuming barycentric coordinates)
        // This would require additional triangle projection logic

        //T f_candidate = sdf.samplePoint(x_candidate);
        T fCandidate = valueAtProjection(*(params.grid), xCandidate,
                                         *(info.B_centerTranslation),
                                         *(info.B_centerRotation));

        // Armijo condition for spatial minimization
        if (fCandidate <= f0 + c * alpha * dfDirection) { return xCandidate; }

        alpha *= rho;
        //Early termination
        if (alpha < 1e-10) break;
    }

    //If Armijo fails, return midpoint
    return x_current + 0.5 * direction;
}

template <typename T>
T FrankWolfeBacktracking(T tstart, T tend, const RigidBodyInfo<T>& initialState/*const EigenVector3<T>& p0,
                const EigenVector3<T>& p1, EigenVector3<T>& p2,
                const grid::Grid<T, T>& grid,*/)
{
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

    //TODO  Compute the barycentric coordinates 𝑢, 𝑣, 𝑤 of the starting iterate.
    //Very importantly: For now we incorrectly assume velocity only for the triangles, not the SDF!

    //Okay this is a very important detail, we have our initial vi is not equation
    // from the paper, it is actually the linear velocity of our triangle rigid body!
    // Thus we must see that vti and vi is different and vti uses getvelocityatpoint!

    EigenVector3<T> p0s = ((initialState.A_p0)).eval();
    EigenVector3<T> p1s = ((initialState.A_p1)).eval();
    EigenVector3<T> p2s = ((initialState.A_p2)).eval();

    EigenVector3<T> vi = *(initialState.A_linearVel);
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

    //if (p0Min > p1Min && p0 > p2Min) { p0Min = }

    DistanceAtTimeParams distanceAtTimeParams{.u = u,
                                              .v = v,
                                              .w = w,
                                              .grid = (initialState.B_sdf),
                                              .p0 = (initialState.A_p0),
                                              .p1 = (initialState.A_p1),
                                              .p2 = (initialState.A_p2)};

    float eps = 1e-7;
    size_t maxIterations = 32u;
    DistanceAtPointParams distanceAtPointParams{.grid = initialState.B_sdf};
    EigenVector3<T> xtip1 = EigenVector3<T>(0, 0, 0);
    EigenVector3<T> xti = EigenVector3<T>(0, 0, 0);

    for (size_t i = 0; i < maxIterations; ++i)
    {
        xti = BarycentricInterpolate(u, v, w, ti, initialState).eval();
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
            T direction = -1.0;
            T initial_step = ti - tstart;
            T gradDotDir = gradPhixti.dot(vti);
            T fCurrent = UnsignedDistanceAtTime<T>(initialState,
                                                   distanceAtTimeParams, ti);
            tip1 = armijoTemporalSearch(
                UnsignedDistanceAtTime<T>, distanceAtTimeParams, initialState,
                tstart, ti, ti, fCurrent, direction * gradDotDir);
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
            T direction = di;
            T initial_step = (di < 0) ? (ti - tstart) : (tend - ti);
            T gradDotDir = gradPhixti.dot(vti) * direction;
            if (di < 0)
            {

                T fCurrent = SignedDistanceAtTime<T>(initialState,
                                                     distanceAtTimeParams, ti);
                tip1 = armijoTemporalSearch(UnsignedDistanceAtTime<T>,
                                            distanceAtTimeParams, initialState,
                                            tstart, ti, ti, fCurrent,
                                            direction * gradDotDir);
            }
            else
            {

                T fCurrent = SignedDistanceAtTime<T>(initialState,
                                                     distanceAtTimeParams, ti);
                tip1 = armijoTemporalSearch(UnsignedDistanceAtTime<T>,
                                            distanceAtTimeParams, initialState,
                                            ti, tend, ti, fCurrent, gradDotDir);
            }
        }
        // Solve spatial sub-problem
        //Note, xtip1 = x_{t_{i+1}}

        xtip1 = BarycentricInterpolate(u, v, w, tip1, initialState).eval();
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
            *(initialState.A_angularVel), p0s, tip1);
        EigenVector3<T> p1_at_ti = getTriangleVertexPosAt(
            *(initialState.A_centerTranslation), *(initialState.A_linearVel),
            *(initialState.A_angularVel), p1s, tip1);
        EigenVector3<T> p2_at_ti = getTriangleVertexPosAt(
            *(initialState.A_centerTranslation), *(initialState.A_linearVel),
            *(initialState.A_angularVel), p2s, tip1);

        T p0Mins = (p0_at_ti).dot(gradPhixtip1);
        T p1Mins = (p1_at_ti).dot(gradPhixtip1);
        T p2Mins = (p2_at_ti).dot(gradPhixtip1);

        if (phixtip1 <= 0) { end = std::min<T>(tip1, end); }
        EigenVector3<T> si;
        //Computing support vertex si
        if (p0Mins <= p1Mins && p0Mins <= p2Mins) { si = p0_at_ti; }
        else if (p1Mins <= p2Mins && p1Mins <= p0Mins) { si = p1_at_ti; }
        else { si = p2_at_ti; }
        xtip1 = armijoSpatialSearch(SignedDistanceAtPoint<T>,
                                    distanceAtPointParams, initialState, xtip1,
                                    si);
        EigenVector3<T> p0_at_tip1 = getTriangleVertexPosAt(
            *(initialState.A_centerTranslation), *(initialState.A_linearVel),
            *(initialState.A_angularVel), p0s, tip1);
        EigenVector3<T> p1_at_tip1 = getTriangleVertexPosAt(
            *(initialState.A_centerTranslation), *(initialState.A_linearVel),
            *(initialState.A_angularVel), p1s, tip1);
        EigenVector3<T> p2_at_tip1 = getTriangleVertexPosAt(
            *(initialState.A_centerTranslation), *(initialState.A_linearVel),
            *(initialState.A_angularVel), p2s, tip1);
        barycentric(p0_at_tip1, p1_at_tip1, p2_at_tip1, xtip1, u, v, w);

        projectToTriangle(u, v, w);
        distanceAtTimeParams.u = u;
        distanceAtTimeParams.v = v;
        distanceAtTimeParams.w = w;

        T phixtip1_2 = valueAtProjection(*(initialState.B_sdf), xtip1,
                                         *(initialState.B_centerTranslation),
                                         *(initialState.B_centerRotation));
        if (std::abs(tip1 - ti) <= eps
            && (std::abs(xtip1.x() - xti.x()) <= eps
                && std::abs(xtip1.y() - xti.y()) <= eps
                && std::abs(xtip1.z() - xti.z()) <= eps)
            && phixtip1_2 >= T(0))
        {
            if (phixtip1_2 >= eps)
            {
                ti = T(0.01);
                tip1 = T(0.01);
            }
            break;
        }
        ti = tip1;
        xti = xtip1;
        /*std::cerr << "ti , tip1 = (" << ti << ", " << tip1 << ")\n";

        std::cerr << "ENDED UP WITH xtip1 = (" << xtip1.x() << ", " << xtip1.y()
                  << ", " << xtip1.z() << ") and xti=" << "(" << xti.x() << ", "
                  << xti.y() << ", " << xti.z() << ")\n";*/
    }
    std::cerr << "ENDED UP WITH tip1 = " << tip1 << " and ti = " << ti << "\n";
    std::cerr << "Ended up with a distance of (from xtip1) "
              << valueAtProjection(*(initialState.B_sdf), xtip1,
                                   *(initialState.B_centerTranslation),
                                   *(initialState.B_centerRotation))
              << " from the solution.\n";
    std::cerr << "Ended up with a distance of (from xti) "
              << valueAtProjection(*(initialState.B_sdf), xti,
                                   *(initialState.B_centerTranslation),
                                   *(initialState.B_centerRotation))
              << " from the solution.\n";

    if (ti <= 0.00500000001)
    {
        //THis code forces debug breakpoint, DELETE LATER when I figure out why TOI=0
        T a = 0;
        T b = a + ti;
        std::cerr << b << ";";
    }
    return std::min<T>(tip1, ti);
}

template <typename T>
T FrankWolfeBacktracking_BENCHMARK_TIME(T tstart, T tend,
                                        const RigidBodyInfo<T>& initialState,
                                        EigenVector3<T>& xtiPoint,
                                        std::vector<T>& minimizerSteps)
{
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

    //TODO  Compute the barycentric coordinates 𝑢, 𝑣, 𝑤 of the starting iterate.
    //Very importantly: For now we incorrectly assume velocity only for the triangles, not the SDF!

    //Okay this is a very important detail, we have our initial vi is not equation
    // from the paper, it is actually the linear velocity of our triangle rigid body!
    // Thus we must see that vti and vi is different and vti uses getvelocityatpoint!

    EigenVector3<T> p0s = ((initialState.A_p0)).eval();
    EigenVector3<T> p1s = ((initialState.A_p1)).eval();
    EigenVector3<T> p2s = ((initialState.A_p2)).eval();

    EigenVector3<T> vi = *(initialState.A_linearVel);
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

    //if (p0Min > p1Min && p0 > p2Min) { p0Min = }

    DistanceAtTimeParams distanceAtTimeParams{.u = u,
                                              .v = v,
                                              .w = w,
                                              .grid = (initialState.B_sdf),
                                              .p0 = (initialState.A_p0),
                                              .p1 = (initialState.A_p1),
                                              .p2 = (initialState.A_p2)};

    float eps = 1e-7;
    size_t maxIterations = 32u;
    DistanceAtPointParams distanceAtPointParams{.grid = initialState.B_sdf};
    EigenVector3<T> xtip1 = EigenVector3<T>(0, 0, 0);
    EigenVector3<T> xti = EigenVector3<T>(0, 0, 0);

    for (size_t i = 0; i < maxIterations; ++i)
    {
        xti = BarycentricInterpolate(u, v, w, ti, initialState).eval();
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
            T direction = -1.0;
            T initial_step = ti - tstart;
            T gradDotDir = gradPhixti.dot(vti);
            T fCurrent = UnsignedDistanceAtTime<T>(initialState,
                                                   distanceAtTimeParams, ti);

            auto gss_start = std::chrono::high_resolution_clock::now();
            tip1 = armijoTemporalSearch(
                UnsignedDistanceAtTime<T>, distanceAtTimeParams, initialState,
                tstart, ti, ti, fCurrent, direction * gradDotDir);
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
            T direction = di;
            T initial_step = (di < 0) ? (ti - tstart) : (end - ti);
            T gradDotDir = gradPhixti.dot(vti) * direction;
            if (di < 0)
            {

                T fCurrent = SignedDistanceAtTime<T>(initialState,
                                                     distanceAtTimeParams, ti);
                auto gss_start = std::chrono::high_resolution_clock::now();
                tip1 = armijoTemporalSearch(UnsignedDistanceAtTime<T>,
                                            distanceAtTimeParams, initialState,
                                            tstart, ti, ti, fCurrent,
                                            direction * gradDotDir);
                auto gss_end = std::chrono::high_resolution_clock::now();
                auto gss_us
                    = std::chrono::duration_cast<std::chrono::nanoseconds>(
                          gss_end - gss_start)
                          .count();
                minimizerSteps.push_back(T(gss_us));
            }
            else
            {

                T fCurrent = SignedDistanceAtTime<T>(initialState,
                                                     distanceAtTimeParams, ti);
                auto gss_start = std::chrono::high_resolution_clock::now();
                tip1 = armijoTemporalSearch(UnsignedDistanceAtTime<T>,
                                            distanceAtTimeParams, initialState,
                                            ti, end, ti, fCurrent, gradDotDir);
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

        xtip1 = BarycentricInterpolate(u, v, w, tip1, initialState).eval();
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
            *(initialState.A_angularVel), p0s, tip1);
        EigenVector3<T> p1_at_ti = getTriangleVertexPosAt(
            *(initialState.A_centerTranslation), *(initialState.A_linearVel),
            *(initialState.A_angularVel), p1s, tip1);
        EigenVector3<T> p2_at_ti = getTriangleVertexPosAt(
            *(initialState.A_centerTranslation), *(initialState.A_linearVel),
            *(initialState.A_angularVel), p2s, tip1);

        T p0Mins = (p0_at_ti).dot(gradPhixtip1);
        T p1Mins = (p1_at_ti).dot(gradPhixtip1);
        T p2Mins = (p2_at_ti).dot(gradPhixtip1);

        if (phixtip1 <= 0) { end = std::min<T>(tip1, end); }
        EigenVector3<T> si;
        //Computing support vertex si
        if (p0Mins <= p1Mins && p0Mins <= p2Mins) { si = p0_at_ti; }
        else if (p1Mins <= p2Mins && p1Mins <= p0Mins) { si = p1_at_ti; }
        else { si = p2_at_ti; }
        auto gss_start = std::chrono::high_resolution_clock::now();
        xtip1 = armijoSpatialSearch(SignedDistanceAtPoint<T>,
                                    distanceAtPointParams, initialState, xtip1,
                                    si);
        auto gss_end = std::chrono::high_resolution_clock::now();
        auto gss_us = std::chrono::duration_cast<std::chrono::nanoseconds>(
                          gss_end - gss_start)
                          .count();
        minimizerSteps.push_back(T(gss_us));
        EigenVector3<T> p0_at_tip1 = getTriangleVertexPosAt(
            *(initialState.A_centerTranslation), *(initialState.A_linearVel),
            *(initialState.A_angularVel), p0s, tip1);
        EigenVector3<T> p1_at_tip1 = getTriangleVertexPosAt(
            *(initialState.A_centerTranslation), *(initialState.A_linearVel),
            *(initialState.A_angularVel), p1s, tip1);
        EigenVector3<T> p2_at_tip1 = getTriangleVertexPosAt(
            *(initialState.A_centerTranslation), *(initialState.A_linearVel),
            *(initialState.A_angularVel), p2s, tip1);
        barycentric(p0_at_tip1, p1_at_tip1, p2_at_tip1, xtip1, u, v, w);

        projectToTriangle(u, v, w);
        distanceAtTimeParams.u = u;
        distanceAtTimeParams.v = v;
        distanceAtTimeParams.w = w;

        T phixtip1_2 = valueAtProjection(*(initialState.B_sdf), xtip1,
                                         *(initialState.B_centerTranslation),
                                         *(initialState.B_centerRotation));
        if (std::abs(tip1 - ti) <= eps
            && (std::abs(xtip1.x() - xti.x()) <= eps
                && std::abs(xtip1.y() - xti.y()) <= eps
                && std::abs(xtip1.z() - xti.z()) <= eps)
            && phixtip1_2 >= T(0))
        {
            if (phixtip1_2 >= eps)
            {
                ti = T(0.01);
                tip1 = T(0.01);
            }
            break;
        }
        ti = tip1;
        xti = xtip1;
    }

    xtiPoint = xtip1;
    return std::min<T>(tip1, ti);
}
} // namespace grid

#endif
