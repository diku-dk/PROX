#ifndef GRID_CCD_FRANK_WOLFE_GOLDEN_SECTION_SEARCH_HPP
#define GRID_CCD_FRANK_WOLFE_GOLDEN_SECTION_SEARCH_HPP

#include <grid_grid.h>
#include <grid_local_optimization.hpp>
#include <chrono>

#include <grid_ccd_frank_wolfe_utils.hpp>

namespace grid
{
template <typename T, typename F>
T GSSMinimize_WHAT_MODIFIED(T lstart, T lend, F func,
                            const DistanceAtTimeParams<T>& params,
                            const RigidBodyInfo<T>& info)
{
    //We init variables
    T phiInv = ((sqrt(T(5)) - T(1)) * T(0.5));
    T r = phiInv;
    T rInv = 1 - r;
    T alpha0 = 0;
    T alpha1 = rInv;
    T alpha2 = r;
    T alpha3 = 1;

    T l0 = lstart;
    T l1 = std::lerp<T>(lstart, lend, alpha1);
    T l2 = std::lerp<T>(lstart, lend, alpha2);
    T l3 = lend;
    T f0 = func(info, params, l0);
    T f1 = func(info, params, l1);
    T f2 = func(info, params, l2);
    T f3 = func(info, params, l3);
    uint16_t it = 0;
    while ((l3 - l0) > 1e-8 /*tol(l1 + l2) && it < 64*/)
    {
        if (std::min<T>(f0, f1) < std::min<T>(f2, f3))
        {
            alpha3 = alpha2;
            l3 = l2;
            f3 = f2;
            alpha2 = alpha1;
            l2 = l1;
            f2 = f1;
            alpha1 = r * alpha2 + rInv * alpha0;
            l1 = std::lerp<T>(lstart, lend, alpha1);
            f1 = func(info, params, l1);
        }
        else
        {
            alpha0 = alpha1;
            l0 = l1;
            f0 = f1;
            alpha1 = alpha2;
            l1 = l2;
            f1 = f2;
            alpha2 = r * alpha1 + rInv * alpha3;
            l2 = std::lerp<T>(lstart, lend, alpha2);
            f2 = func(info, params, l2);
        }
        it++;
    }
    T lmid = std::lerp<T>(lstart, lend, 0.5 * (alpha0 + alpha3));
    T fmid = func(info, params, lmid);
    T lmin;
    if (f0 < fmid && f0 < f3) { lmin = l0; }
    else if (fmid < f3) { lmin = lmid; }
    else { lmin = l3; }
    return lmin;
}

/*template <typename T, typename F>
EigenVector3<T> GSSMinimize_WHAT_MODIFIED(
    T lstart, T lend, EigenVector3<T> pstart, EigenVector3<T> pend, F func,
    const DistanceAtPointParams<T>& params, const RigidBodyInfo<T>& info)
{
    T diff = (lstart - lend);
    if (std::abs<T>(diff) < 1e-8) { return pstart; }
    //We init variables
    T phiInv = ((sqrt(T(5)) - T(1)) * T(0.5));
    T r = phiInv;
    T rInv = 1 - r;
    T alpha0 = 0;
    T alpha1 = rInv;
    T alpha2 = r;
    T alpha3 = 1;

    T l0 = lstart;
    T l1 = std::lerp<T>(lstart, lend, alpha1);
    T l2 = std::lerp<T>(lstart, lend, alpha2);
    T l3 = lend;
    T f0 = func(info, params, l0, pstart, pend);
    T f1 = func(info, params, l1, pstart, pend);
    T f2 = func(info, params, l2, pstart, pend);
    T f3 = func(info, params, l3, pstart, pend);
    uint16_t it = 0;
    while ((l3 - l0) <= tol(l1 + l2) && it < 8)
    {
        if (std::min<T>(f0, f1) < std::min<T>(f2, f3))
        {
            alpha3 = alpha2;
            l3 = l2;
            f3 = f2;
            alpha2 = alpha1;
            l2 = l1;
            f2 = f1;
            alpha1 = r * alpha2 + rInv * alpha0;
            l1 = std::lerp<T>(lstart, lend, alpha1);
            f1 = func(info, params, l1, pstart, pend);
        }
        else
        {
            alpha0 = alpha1;
            l0 = l1;
            f0 = f1;
            alpha1 = alpha2;
            l1 = l2;
            f1 = f2;
            alpha2 = r * alpha1 + rInv * alpha3;
            l2 = std::lerp<T>(lstart, lend, alpha2);
            f2 = func(info, params, l2, pstart, pend);
        }
        it++;
    }
    T lmid = std::lerp<T>(lstart, lend, 0.5 * (alpha0 + alpha3));
    T fmid = func(info, params, lmid, pstart, pend);
    T lmin;
    if (f0 < fmid && f0 < f3) { lmin = l0; }
    else if (fmid < f3) { lmin = lmid; }
    else { lmin = l3; }
    return pstart + (pend - pstart) * lmin;
}*/

/*template <typename T, typename F>
T GSSMinimize_WHAT(T lstart, T lend, F func,
                   const DistanceAtTimeParams<T>& params,
                   const RigidBodyInfo<T>& info)
{
    //We init variables
    T phiInv = ((sqrt(T(5)) - T(1)) * T(0.5));
    T r = phiInv;
    T rInv = 1.0 - r;
    T alpha0 = 0;
    T alpha1 = rInv;
    T alpha2 = r;
    T alpha3 = 1;

    T l0 = lstart;
    T l1 = std::lerp<T>(lstart, lend, alpha1);
    T l2 = std::lerp<T>(lstart, lend, alpha2);
    T l3 = lend;
    T f0 = func(info, params, l0);
    T f1 = func(info, params, l1);
    T f2 = func(info, params, l2);
    T f3 = func(info, params, l3);
    uint16_t it = 0;
    while ((l3 - l0) <= tol(l1 + l2) && it < 32)
    {
        if (std::min<T>(f0, f1) < std::min<T>(f2, f3))
        {
            alpha3 = alpha2;
            l3 = l2;
            f3 = f2;
            alpha2 = alpha1;
            l2 = l1;
            f2 = f1;
            alpha1 = r * alpha2 + rInv * alpha0;
            l1 = std::lerp<T>(lstart, lend, alpha1);
            f1 = func(info, params, l1);
        }
        else
        {
            alpha0 = alpha1;
            l0 = l1;
            f0 = f1;
            alpha1 = alpha2;
            l1 = l2;
            f1 = f2;
            alpha2 = r * alpha1 + rInv * alpha3;
            l2 = std::lerp<T>(lstart, lend, alpha2);
            f2 = func(info, params, l2);
        }
        it++;
    }
    T lmid = std::lerp<T>(lstart, lend, 0.5 * (alpha0 + alpha3));
    T fmid = func(info, params, lmid);
    T lmin;
    if (f0 < fmid && f0 < f3) { lmin = l0; }
    else if (fmid < f3) { lmin = lmid; }
    else { lmin = l3; }
    return lmin;
}*/

template <typename T, typename F>
EigenVector3<T> GSSMinimize_WHAT(EigenVector3<T> lstart, EigenVector3<T> lend,
                                 F func, const DistanceAtPointParams<T>& params,
                                 const RigidBodyInfo<T>& info)
{
    //We init variables
    T phiInv = ((sqrt(T(5)) - T(1)) * T(0.5));
    T r = phiInv;
    T rInv = 1 - r;
    T alpha0 = 0;
    T alpha1 = rInv;
    T alpha2 = r;
    T alpha3 = 1;

    EigenVector3<T> l0 = lstart;
    EigenVector3<T> l1 = lerp<T>(lstart, lend, alpha1);
    EigenVector3<T> l2 = lerp<T>(lstart, lend, alpha2);
    EigenVector3<T> l3 = lend;
    T f0 = func(info, params, l0);
    T f1 = func(info, params, l1);
    T f2 = func(info, params, l2);
    T f3 = func(info, params, l3);
    T tol = 1e-5;
    //TODO NOT IMPLEMENTED CORRECTLY
    //while ((l3 - l0) <= tol((l1 + l2)))
    uint16_t it = 0;
    while ((l3 - l0).norm() > tol /*((l1 + l2).norm() * tol)*/ /* && it < 8*/)
    {
        if (std::min<T>(f0, f1) < std::min<T>(f2, f3))
        {
            alpha3 = alpha2;
            l3 = l2;
            f3 = f2;
            alpha2 = alpha1;
            l2 = l1;
            f2 = f1;
            alpha1 = r * alpha2 + rInv * alpha0;
            l1 = lerp<T>(lstart, lend, alpha1);
            f1 = func(info, params, l1);
        }
        else
        {
            alpha0 = alpha1;
            l0 = l1;
            f0 = f1;
            alpha1 = alpha2;
            l1 = l2;
            f1 = f2;
            alpha2 = r * alpha1 + rInv * alpha3;
            l2 = lerp<T>(lstart, lend, alpha2);
            f2 = func(info, params, l2);
        }
        ++it;
    }
    EigenVector3<T> lmid = lerp<T>(lstart, lend, T(0.5 * (alpha0 + alpha3)));
    T fmid = func(info, params, lmid);
    EigenVector3<T> lmin;
    if (f0 < fmid && f0 < f3) { lmin = l0; }
    else if (fmid < f3) { lmin = lmid; }
    else { lmin = l3; }
    return lmin;
}

template <typename T, typename F>
T GSSMinimize(T a, T b, F func, const DistanceAtTimeParams<T>& params,
              const RigidBodyInfo<T>& info, T tolerance = T(1e-5))
{
    if (std::abs(a - b) < 1e-8) { return a; }
    const T phi = (1 + std::sqrt(5)) / 2;
    const T invphi = 1 / phi;

    // initial interior points
    T c = b - (b - a) / phi;
    T d = a + (b - a) / phi;

    T f_c = func(info, params, c);
    T f_d = func(info, params, d);
    int it = 0;

    while ((b - a) > tolerance /* && it < 8*/)
    {
        if (f_c < f_d)
        {
            b = d;
            d = c;
            f_d = f_c;
            c = b - (b - a) / phi;
            f_c = func(info, params, c);
        }
        else
        {
            a = c;
            c = d;
            f_c = f_d;
            d = a + (b - a) / phi;
            f_d = func(info, params, d);
        }
        it++;
    }

    //Return midpoint!
    return 0.5 * (a + b);
}

/*template <typename T, typename F>
EigenVector3<T> GSSMinimize(EigenVector3<T> a, EigenVector3<T> b, F func,
                            const DistanceAtPointParams<T>& params,
                            const RigidBodyInfo<T>& info, T tolerance = T(1e-5))
{
    const T phi = (1 + std::sqrt(5)) / 2;
    const T invphi = 1 / phi;

    // initial interior points
    EigenVector3<T> c = b - (b - a) / phi;
    EigenVector3<T> d = a + (b - a) / phi;

    T f_c = func(info, params, c);
    T f_d = func(info, params, d);
    int it = 0;

    while (((b - a).x() > tolerance && (b - a).y() > tolerance
            && (b - a).z() > tolerance)
           && it < 8)
    {
        if (f_c < f_d)
        {
            b = d;
            d = c;
            f_d = f_c;
            c = b - (b - a) / phi;
            f_c = func(info, params, c);
        }
        else
        {
            a = c;
            c = d;
            f_c = f_d;
            d = a + (b - a) / phi;
            f_d = func(info, params, d);
        }
        it++;
    }

    //Return midpoint!
    return 0.5 * (a + b);
}*/

template <typename T>
T FrankWolfeGSS(T tstart, T tend, const RigidBodyInfo<T>& initialState/*const EigenVector3<T>& p0,
                const EigenVector3<T>& p1, EigenVector3<T>& p2,
                const grid::Grid<T, T>& grid,*/)
{
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
        EigenVector3<T> gradPhixti
            = gradientAtProjection(xti, *(initialState.B_sdf),
                                   *(initialState.B_centerTranslation),
                                   *(initialState.B_centerRotation))
                  .normalized();
        if (phixti <= 0)
        {
            end = std::min<T>(ti, end);
            tip1 = GSSMinimize_WHAT_MODIFIED(
                tstart, ti, UnsignedDistanceAtTime<T>, distanceAtTimeParams,
                initialState);
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
                tip1 = GSSMinimize_WHAT_MODIFIED(
                    tstart, ti, SignedDistanceAtTime<T>, distanceAtTimeParams,
                    initialState);
            }
            else
            {
                tip1 = GSSMinimize_WHAT_MODIFIED(
                    ti, end, SignedDistanceAtTime<T>, distanceAtTimeParams,
                    initialState);
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
        /*        T p0Min = (*(initialState.A_p0)).dot(gradPhixtip1);
        T p1Min = (*(initialState.A_p1)).dot(gradPhixtip1);
        T p2Min = (*(initialState.A_p2)).dot(gradPhixtip1);

        if (phixtip1 <= 0) { end = std::min<T>(tip1, end); }
        EigenVector3<T> si;
        //Computing support vertex pi
        if (p0Min <= p1Min && p0Min <= p2Min) { si = *(initialState.A_p0); }
        else if (p1Min <= p2Min && p1Min <= p0Min)
        {
            si = *(initialState.A_p1);
        }
        else { si = *(initialState.A_p2); }*/
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
        xtip1 = GSSMinimize_WHAT(xtip1, si, SignedDistanceAtPoint<T>,
                                 distanceAtPointParams, initialState);
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

        //projectToTriangle(u, v, w);
        distanceAtTimeParams.u = u;
        distanceAtTimeParams.v = v;
        distanceAtTimeParams.w = w;
        /*if (u > 1.01 || v > 1.01 || w > 1.01)
        {
            std::cerr << u << ", " << v << ", " << w << "\n";
            throw std::runtime_error("TOO HIGH");
        }
        if (u < -0.01 || v < -0.01 || w < -0.01)
        {
            std::cerr << u << ", " << v << ", " << w << "\n";
            throw std::runtime_error("TOO LOW");
        }*/

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
    /*std::cerr << "ENDED UP WITH tip1 = " << tip1 << " and ti = " << ti << "\n";
    std::cerr << "Ended up with a distance of (from xtip1) "
              << valueAtProjection(*(initialState.B_sdf), xtip1,
                                   *(initialState.B_centerTranslation),
                                   *(initialState.B_centerRotation))
              << " from the solution.\n";
    std::cerr << "Ended up with a distance of (from xti) "
              << valueAtProjection(*(initialState.B_sdf), xti,
                                   *(initialState.B_centerTranslation),
                                   *(initialState.B_centerRotation))
              << " from the solution.\n";*/

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
T FrankWolfeGSS_BENCHMARK_TIME(T tstart, T tend,
                               const RigidBodyInfo<T>& initialState,
                               EigenVector3<T>& xtiPoint,
                               std::vector<T>& minimizerSteps)
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
        xti = BarycentricInterpolate(distanceAtTimeParams.u,
                                     distanceAtTimeParams.v,
                                     distanceAtTimeParams.w, ti, initialState)
                  .eval();
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

        xtip1 = BarycentricInterpolate(
                    distanceAtTimeParams.u, distanceAtTimeParams.v,
                    distanceAtTimeParams.w, tip1, initialState)
                    .eval();
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
        barycentric(p0_at_tip1, p1_at_tip1, p2_at_tip1, xtip1,
                    distanceAtTimeParams.u, distanceAtTimeParams.v,
                    distanceAtTimeParams.w);

        projectToTriangle(distanceAtTimeParams.u, distanceAtTimeParams.v,
                          distanceAtTimeParams.w);
        /*distanceAtTimeParams.u = u;
        distanceAtTimeParams.v = v;
        distanceAtTimeParams.w = w;*/

        T phixtip1_2 = valueAtProjection(*(initialState.B_sdf), xtip1,
                                         *(initialState.B_centerTranslation),
                                         *(initialState.B_centerRotation));
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

    xtiPoint = xtip1;
    return std::min<T>(tip1, ti);
}
}

#endif
