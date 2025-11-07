#ifndef GRID_CCD_FRANK_WOLFE_CONVERGENCE_PLOTS_HPP
#define GRID_CCD_FRANK_WOLFE_CONVERGENCE_PLOTS_HPP

#include <grid_ccd_frank_wolfe_include_all.hpp>

namespace grid
{
template <typename T>
T FrankWolfeBRENT_BENCHMARK_TIME_CONVERGENCE(
    T tstart, T tend, const RigidBodyInfo<T>& initialState,
    EigenVector3<T>& xtiPoint, std::vector<T>& minimizerSteps,
    size_t& iterations, std::vector<T>& xtixiErr, std::vector<T>& tip1tiErr)
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
        iterations = 1;
        xtixiErr.push_back(T(0));
        tip1tiErr.push_back(T(0));
        return tend;
    }
    if (tstart >= tend)
    {
        iterations = 1;
        xtixiErr.push_back(T(0));
        tip1tiErr.push_back(T(0));
        return tend;
    }

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

    float eps = 1e-6;
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
        EigenVector3<T> diff = (xtip1 - xti);
        T onextixiErr = std::abs<T>(diff.norm());
        T onetip1tiErr = std::abs<T>(tip1 - ti);
        xtixiErr.push_back(onextixiErr);
        tip1tiErr.push_back(onetip1tiErr);
        if (std::abs(tip1 - ti) <= eps
            && (std::abs(xtip1.x() - xti.x()) <= eps
                && std::abs(xtip1.y() - xti.y()) <= eps
                && std::abs(xtip1.z() - xti.z()) <= eps)
            && phixtip1_2 >= -eps)
        {
            /*if (phixtip1_2 >= eps)
            {
                ti = T(0.01);
                tip1 = T(0.01);
            }*/
            break;
        }
        if (its > maxIterations && phixtip1_2 >= -eps) { break; }
        /*if (its > maxIterations)
        {
            std::cerr << "CUrrent it: " << i << " and xtip-value: ("
                      << xtip1.x() << ", " << xtip1.y() << ", " << xtip1.z()
                      << ") and phixtip1_2 =" << phixtip1_2 << "\n";
        }*/
        ti = tip1;
        xti = xtip1;
    }

    iterations = its;

    xtiPoint = xtip1;
    if (tip1tiErr.size() != xtixiErr.size())
    {
        throw std::runtime_error("DIFF!");
    }
    return std::min<T>(tip1, ti);
}
} // namespace grid

namespace grid
{
template <typename T>
T FrankWolfeGSS_BENCHMARK_TIME_CONVERGENCE(
    T tstart, T tend, const RigidBodyInfo<T>& initialState,
    EigenVector3<T>& xtiPoint, std::vector<T>& minimizerSteps,
    size_t& iterations, std::vector<T>& xtixiErr, std::vector<T>& tip1tiErr)
{
    //return FrankWolfeGSS(tstart, tend, initialState);

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
        iterations = 1;
        xtixiErr.push_back(T(0));
        tip1tiErr.push_back(T(0));
        return tend;
    }

    EigenVector3<T> p0s = ((initialState.A_p0)).eval();
    EigenVector3<T> p1s = ((initialState.A_p1)).eval();
    EigenVector3<T> p2s = ((initialState.A_p2)).eval();

    EigenVector3<T> vi = *(initialState.A_linearVel);
    EigenVector3<T> gradP0
        = gradientAtProjection(p0s, *(initialState.B_sdf),
                               *(initialState.B_centerTranslation),
                               *(initialState.B_centerRotation))
              .normalized();
    EigenVector3<T> gradP1
        = gradientAtProjection(p1s, *(initialState.B_sdf),
                               *(initialState.B_centerTranslation),
                               *(initialState.B_centerRotation))
              .normalized();
    EigenVector3<T> gradP2
        = gradientAtProjection(p2s, *(initialState.B_sdf),
                               *(initialState.B_centerTranslation),
                               *(initialState.B_centerRotation))
              .normalized();
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

    float eps = 1e-7;
    size_t maxIterations = 32u;
    //The higher the below iterations, the higher likelihood we will get no false negatives
    size_t hardStopMaxIterations = 1000u;
    size_t its = 0;
    DistanceAtPointParams distanceAtPointParams{.grid = initialState.B_sdf};
    EigenVector3<T> xtip1 = EigenVector3<T>(0, 0, 0);
    EigenVector3<T> xti = EigenVector3<T>(0, 0, 0);

    bool standStill = false;
    for (size_t i = 0; i < hardStopMaxIterations; ++i)
    {
        its += 1;
        xti = BarycentricInterpolate(u, v, w, ti, initialState).eval();
        //I assume (but only assumption that we call with ti and xti!
        EigenVector3<T> vti = getVelocityAtPoint(initialState, xti, ti);
        T phixti = valueAtProjection(*(initialState.B_sdf), xti,
                                     *(initialState.B_centerTranslation),
                                     *(initialState.B_centerRotation));
        EigenVector3<T> gradPhixti
            = gradientAtProjection(xti, *(initialState.B_sdf),
                                   *(initialState.B_centerTranslation),
                                   *(initialState.B_centerRotation))
                  .normalized();
        if (phixti <= 0)
        {
            end = std::min<T>(ti, end);
            auto gss_start = std::chrono::high_resolution_clock::now();
            tip1 = GSSMinimize_WHAT_MODIFIED(
                tstart, end, UnsignedDistanceAtTime<T>, distanceAtTimeParams,
                initialState);
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
                tip1 = GSSMinimize_WHAT_MODIFIED(
                    tstart, ti, SignedDistanceAtTime<T>, distanceAtTimeParams,
                    initialState);
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
                tip1 = GSSMinimize_WHAT_MODIFIED(
                    ti, end, SignedDistanceAtTime<T>, distanceAtTimeParams,
                    initialState);
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
        EigenVector3<T> gradPhixtip1
            = gradientAtProjection(xtip1, *(initialState.B_sdf),
                                   *(initialState.B_centerTranslation),
                                   *(initialState.B_centerRotation))
                  .normalized();

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
        T onextixiErr = std::abs<T>((xtip1 - xti).norm());
        T onetip1tiErr = std::abs<T>(tip1 - ti);
        xtixiErr.push_back(onextixiErr);
        tip1tiErr.push_back(onetip1tiErr);
        if (std::abs(tip1 - ti) <= eps
            && (std::abs(xtip1.x() - xti.x()) <= eps
                && std::abs(xtip1.y() - xti.y()) <= eps
                && std::abs(xtip1.z() - xti.z()) <= eps))
        {
            if (phixtip1_2 >= -eps) { break; }
            else
            {
                if (standStill) { break; }
                standStill = true;
            }
        }
        if (its > maxIterations && phixtip1_2 >= -eps) { break; }
        ti = tip1;
        xti = xtip1;
    }
    iterations = its;

    xtiPoint = xtip1;
    return std::min<T>(tip1, ti);
}
} // namespace grid

namespace grid
{
template <typename T>
T FrankWolfeBacktracking_BENCHMARK_TIME_CONVERGENCE(
    T tstart, T tend, const RigidBodyInfo<T>& initialState,
    EigenVector3<T>& xtiPoint, std::vector<T>& minimizerSteps,
    size_t& iterations, std::vector<T>& xtixiErr, std::vector<T>& tip1tiErr)
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
        iterations = 1;
        xtixiErr.push_back(T(0));
        tip1tiErr.push_back(T(0));
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
        T onextixiErr = std::abs<T>((xtip1 - xti).norm());
        T onetip1tiErr = std::abs<T>(tip1 - ti);
        xtixiErr.push_back(onextixiErr);
        tip1tiErr.push_back(onetip1tiErr);
        iterations = i;
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

namespace grid
{
template <typename T>
T performProjectedGradientDescent_CONVERGENCE(
    T tstart, T tend, const RigidBodyInfo<T>& initialState,
    EigenVector3<T>& contactPoint, std::vector<T>& minimizerSteps,
    size_t& iterations, std::vector<T>& xtixiErr, std::vector<T>& tip1tiErr)
{
    T minVel = 0.0000001;
    if ((*(initialState.A_linearVel)).norm() < minVel
        && (*(initialState.A_angularVel)).norm() < minVel)
    {
        iterations = 1;
        xtixiErr.push_back(T(0));
        tip1tiErr.push_back(T(0));
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
    uint32_t maxIts = 1024;
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
        }
        //alpha = (T(1.0) / (std::abs<T>(dt) / tend)) * 0.0001;
        //alpha = 1e-6;

        EigenVector4<T> tmpSolution = initialSolution - alpha * gradient;
        ProjectToSimplex(tmpSolution.x(), tmpSolution.y(), tmpSolution.z());
        tmpSolution.w() = std::clamp<T>(tmpSolution.w(), tstart, tend);
        initialSolution = tmpSolution;

        T onextixiErr
            = std::abs<T>((initialSolution - initialSolutionLast).norm());
        T onetip1tiErr = onextixiErr;
        xtixiErr.push_back(onextixiErr);
        tip1tiErr.push_back(onetip1tiErr);
        iterations = iteration;

        if ((initialSolution - initialSolutionLast).norm() < 1e-7) { break; }
        initialSolutionLast = initialSolution;
    }
    if (std::abs<T>(phival) < tol) { ; }
    contactPoint = BarycentricInterpolate(
        initialSolution.x(), initialSolution.y(), initialSolution.z(),
        initialSolution.w(), initialState);
    return initialSolution.w();

    //    return tend; // No collision found
}
} //namespace grid

#endif // GRID_CCD_FRANK_WOLFE_CONVERGENCE_PLOTS_HPP
