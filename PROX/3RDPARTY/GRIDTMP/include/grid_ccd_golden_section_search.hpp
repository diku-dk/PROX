#ifndef GRID_CCD_GOLDEN_SECTION_SEARCH_HPP
#define GRID_CCD_GOLDEN_SECTION_SEARCH_HPP

#include <grid_grid.h>
#include <grid_local_optimization.hpp>

namespace grid
{
/*enum FuncType
{
    SignedDistanceAtTime,
    UnsignedDistanceAtTime,
    SignedDistanceAtPoint
};

template <typename T> T someFunction(T val) { return val + val; }*/

template <typename T> T tol(T val)
{
    T tol = 1e-5;
    //TODO NOT IMPLEMENTED CORRECTLY
    return tol * val;
}

template <typename T> struct RigidBodyInfo
{
    //For rigid body A:
    const EigenVector3<T>* A_p0;
    const EigenVector3<T>* A_p1;
    const EigenVector3<T>* A_p2;
    const EigenVector3<T>* A_linearVel;
    const EigenVector3<T>* A_angularVel;
    const EigenVector3<T>* A_centerTranslation;

    const grid::Grid<T, T>* B_sdf;
    const EigenVector3<T>* B_linearVel;
    const EigenVector3<T>* B_angularVel;
    const EigenVector3<T>* B_centerTranslation;
};

template <typename T> struct TriangleAtTimeInfo
{
    const EigenVector3<T> A_p0;
    const EigenVector3<T> A_p1;
    const EigenVector3<T> A_p2;
};

template <typename T>
TriangleAtTimeInfo<T> getTriangleAtTime(T t,
                                        const RigidBodyInfo<T>& initialState)
{
    //Current center position due to linear motion
    EigenVector3<T> currentCenter
        = initialState.position + initialState.linearVelocity * t;

    T angle = initialState.angularVelocity.norm() * t;

    EigenQuaternion<T> rotation;
    if (angle > 1e-10)
    {
        EigenVector3<T> axis = initialState.angularVelocity.normalized();
        rotation = EigenQuaternion<T>(Eigen::AngleAxis<T>(angle, axis));
    }
    else { rotation = EigenQuaternion<T>::Identity(); }

    // Rotate and translate each vertex
    TriangleAtTimeInfo tInfo{
        .A_p0 = currentCenter + rotation.rotateVector(initialState.a_local),
        .A_p1 = currentCenter + rotation.rotateVector(initialState.b_local),
        .A_p2 = currentCenter + rotation.rotateVector(initialState.c_local)};

    return tInfo;
}

template <typename T>
EigenVector3<T> getVelocityAtPoint(const RigidBodyInfo<T>& initialState,
                                   const EigenVector3<T>& point, T t)
{
    //Get the center g of rotation at time t
    EigenVector3<T> currentCenter
        = initialState.position + initialState.linearVelocity * t;

    //Paper's formula: v_t = v_g + ω_g × (x_t - g)
    EigenVector3<T> radius = point - currentCenter;
    EigenVector3<T> velocity = initialState.linearVelocity
                             + initialState.angularVelocity.cross(radius);

    return velocity;
}

//I understand the ti componenet as the triangle at time ti with the u,v,w interpolation
/*template <typename T>
EigenVector3<T>
BarycentricInterpolate(T u, T v, T w, T ti, const EigenVector3<T>& p0,
                       const EigenVector3<T>& p1, const EigenVector3<T>& p2)
{
    return u * p0 + v * p1 + v * p2;
}*/
template <typename T>
EigenVector3<T> BarycentricInterpolate(T u, T v, T w, T ti,
                                       const RigidBodyInfo<T>& info)
{
    TriangleAtTimeInfo<T> triangle = getTriangleAtTime(ti, info);
    return triangle.A_p0 * u + triangle.A_p1 * v + triangle.A_p2 * w;
}

template <typename T> struct DistanceAtTimeParams
{
    const T u;
    const T v;
    const T w;
    const grid::Grid<T, T>* grid;
    const EigenVector3<T>* p0;
    const EigenVector3<T>* p1;
    const EigenVector3<T>* p2;
};

template <typename T>
T UnsignedDistanceAtTime(const RigidBodyInfo<T>& info,
                         const DistanceAtTimeParams<T>& params, T t)
{
    EigenVector3<T> x
        = BarycentricInterpolate(params.u, params.v, params.w, t, info);
    return std::abs<T>(grid::value_at_2(params.grid, x));
}

template <typename T>
T SignedDistanceAtTime(const RigidBodyInfo<T>& info,
                       const DistanceAtTimeParams<T>& params, T t)
{
    EigenVector3<T> x
        = BarycentricInterpolate(params.u, params.v, params.w, t, info);
    return (grid::value_at_2(params.grid, x));
}

template <typename T> struct DistanceAtPointParams
{
    const grid::Grid<T, T>* grid;
};

template <typename T>
T SignedDistanceAtPoint(const DistanceAtPointParams<T>& params,
                        EigenVector3<T> x)
{
    return (grid::value_at_2(params.grid, x));
}

template <typename T, typename F>
T GSSMinimize(T lstart, T lend, F func, const DistanceAtTimeParams<T>& params,
              const RigidBodyInfo<T>& info)
{
    //We init variables
    T phiInv = T((sqrt(5) - 1) * 0.5);
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
    while ((l3 - l0) <= tol(l1 + l2))
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
            f2 = func(info, params, l2);
        }
    }
    T lmid = std::lerp<T>(lstart, lend, 0.5 * (alpha0 + alpha3));
    T fmid = func(info, params, lmid);
    T lmin;
    if (f0 < fmid && f0 < f3) { lmin = l0; }
    else if (fmid < f3) { lmin = lmid; }
    else { lmin = l3; }
    return lmin;
}

template <typename T>
EigenVector3<T> lerp(EigenVector3<T> x0, EigenVector3<T> x1, T alpha)
{
    EigenVector3<T> val = EigenVector3<T>(std::lerp<T>(x0.x(), x1.x(), alpha),
                                          std::lerp<T>(x0.y(), x1.y(), alpha),
                                          std::lerp<T>(x0.z(), x1.z(), alpha));
    return val;
}

template <typename T, typename F>
EigenVector3<T> GSSMinimize(EigenVector3<T> lstart, EigenVector3<T> lend,
                            F func, const DistanceAtPointParams<T>& params)
{
    //We init variables
    T phiInv = T((sqrt(5) - 1) * 0.5);
    T r = phiInv;
    T rInv = 1 - r;
    T alpha0 = 0;
    T alpha1 = rInv;
    T alpha2 = r;
    T alpha3 = 1;

    EigenVector3<T> l0 = lstart;
    EigenVector3<T> l1 = lerp(lstart, lend, alpha1);
    EigenVector3<T> l2 = lerp(lstart, lend, alpha2);
    EigenVector3<T> l3 = lend;
    T f0 = func(params, l0);
    T f1 = func(params, l1);
    T f2 = func(params, l2);
    T f3 = func(params, l3);
    while ((l3 - l0) <= tol(l1 + l2))
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
            l1 = lerp(lstart, lend, alpha1);
            f1 = func(params, l1);
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
            f2 = func(params, l2);
        }
    }
    EigenVector3<T> lmid = lerp(lstart, lend, 0.5 * (alpha0 + alpha3));
    T fmid = func(params, lmid);
    EigenVector3<T> lmin;
    if (f0 < fmid && f0 < f3) { lmin = l0; }
    else if (fmid < f3) { lmin = lmid; }
    else { lmin = l3; }
    return lmin;
}

template <typename T>
void computeBarycentricCoordinates(const EigenVector3<T>& p0,
                                   const EigenVector3<T>& p1,
                                   EigenVector3<T>& p2,
                                   const EigenVector3<T>& samplePoint, T& u,
                                   T& v, T& w)
{
    const T x = samplePoint.x();
    const T y = samplePoint.y();

    const T x1 = p0.x();
    const T y1 = p0.y();

    const T x2 = p1.x();
    const T y2 = p1.y();

    const T x3 = p2.x();
    const T y3 = p2.y();

    const T y2y3 = y2 - y3;
    const T x1x3 = x1 - x3;
    const T x3x2 = x3 - x2;
    const T y1y3 = y1 - y3;
    const T y3y1 = y3 - y1;
    const T xx3 = x - x3;
    const T yy3 = y - y3;

    const T detT = y2y3 * x1x3 + x3x2 * y1y3;
    u = (y2y3 * xx3 + x3x2 * yy3) / (detT);
    v = (y3y1 * xx3 + x1x3 * yy3) / (detT);
    w = 1 - u - v;
}

template <typename T>
T FrankWolfeGSS(T tstart, T tend, const RigidBodyInfo<T>& initialState/*const EigenVector3<T>& p0,
                const EigenVector3<T>& p1, EigenVector3<T>& p2,
                const grid::Grid<T, T>& grid,*/)
{
    T t1 = tstart;
    T ti = t1;
    T tip1 = ti;
    //TODO  Compute the barycentric coordinates 𝑢, 𝑣, 𝑤 of the starting iterate.
    //Very importantly: For now we incorrectly assume velocity only for the triangles, not the SDF!
    EigenVector3<T> vi = 0; //TODO VELOCITY
    EigenVector3<T> gradP0
        = grid::computeGradient_Working(initialState.A_p0, initialState.B_sdf);
    EigenVector3<T> gradP1
        = grid::computeGradient_Working(initialState.A_p1, initialState.B_sdf);
    EigenVector3<T> gradP2
        = grid::computeGradient_Working(initialState.A_p2, initialState.B_sdf);
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
                                              .grid = &initialState.B_sdf,
                                              .p0 = &initialState.A_p0,
                                              .p1 = &initialState.A_p1,
                                              .p2 = &initialState.A_p2};

    float eps = 1e-6;
    DistanceAtPointParams distanceAtPointParams{.grid = &initialState.B_sdf};
    size_t maxIterations = 8u;
    for (size_t i = 0; i < maxIterations; ++i)
    {
        EigenVector3<T> xti = BarycentricInterpolate(u, v, w, ti, initialState);
        //I assume (but only assumption that we call with ti and xti!
        vi = getVelocityAtPoint(initialState, xti, ti);
        T phixti = grid::value_at_2(initialState.B_sdf, xti);
        EigenVector3<T> gradPhixti
            = grid::computeGradient_Working(xti, initialState.B_sdf);
        if (phixti <= 0)
        {
            tend = std::min<T>(ti, tend);
            tip1 = GSSMinimize(tstart, ti, UnsignedDistanceAtTime,
                               distanceAtTimeParams, initialState);
        }
        else
        {
            //Compute Direction
            T di = T(-1);
            //Below can be std::copysign not sure
            //            if (phixti > 0) { di = -sign(gradPhixti.dot(vi)); }
            if (phixti > 0) { di = -std::copysign(gradPhixti.dot(vi)); }
            //Direction sign test
            if (di < 0)
            {
                tip1 = GSSMinimize(tstart, ti, SignedDistanceAtTime,
                                   distanceAtTimeParams, initialState);
            }
            else
            {
                tip1 = GSSMinimize(ti, tend, SignedDistanceAtTime,
                                   distanceAtTimeParams, initialState);
            }
        }
        // Solve spatial sub-problem
        EigenVector3<T> xtip1
            = BarycentricInterpolate(u, v, w, tip1, initialState);
        //EigenVector3<T> vtip1 = 0; //TODO compute v_{t_{i+1}}
        EigenVector3<T> vtip1 = getVelocityAtPoint(initialState, xtip1, tip1);
        T phixtip1 = grid::value_at_2(initialState.B_sdf, xtip1);
        EigenVector3<T> gradPhixtip1
            = grid::computeGradient_Working(xtip1, initialState.B_sdf);
        T p0Min = initialState.A_p0.dot(gradPhixtip1);
        T p1Min = initialState.A_p1.dot(gradPhixtip1);
        T p2Min = initialState.A_p2.dot(gradPhixtip1);
        EigenVector3<T> si;
        //Computing support vertex pi
        if (p0Min > p1Min && p0Min > p2Min) { si = initialState.A_p0; }
        else if (p1Min > p2Min && p1Min > p0Min) { si = initialState.A_p1; }
        else { si = initialState.A_p2; }
        xtip1 = GSSMinimize(xtip1, si, SignedDistanceAtPoint,
                            distanceAtPointParams);
        //TODO: Update barycentric coordinates 𝑢, 𝑣, 𝑤 using x®𝑡𝑖+1
        computeBarycentricCoordinates(initialState.A_p0, initialState.A_p1,
                                      initialState.A_p2, xtip1, u, v, w);

        if (std::abs(tip1 - ti) <= eps
            && (std::abs(xtip1.x() - xti.x()) <= eps
                && std::abs(xtip1.y() - xti.y()) <= eps
                && std::abs(xtip1.z() - xti.z()) <= eps))
        {
            break;
        }
    }
}
} // namespace grid
#endif // GRID_CCD_GOLDEN_SECTION_SEARCH_HPP
