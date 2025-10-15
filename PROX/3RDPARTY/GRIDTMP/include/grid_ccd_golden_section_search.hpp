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

//Note in the below we assume the sample point is in world coordinates, but the sdf is in local coords!
template <typename T>
T valueAtProjection(const grid::Grid<T, T>& sdf,
                    const EigenVector3<T>& samplePoint,
                    const EigenVector3<T>& translationSDF,
                    const EigenQuaternion<T>& rotationSDF)
{
    EigenVector3<T> localSamplePoint
        = rotationSDF.inverse() * (samplePoint - translationSDF);
    return grid::value_at_2(sdf, localSamplePoint);
}

template <typename T>
EigenVector3<T> gradientAtProjection(const EigenVector3<T>& samplePoint,
                                     const grid::Grid<T, T>& sdf,
                                     const EigenVector3<T>& translationSDF,
                                     const EigenQuaternion<T>& rotationSDF)
{
    EigenVector3<T> localSamplePoint
        = rotationSDF.inverse() * (samplePoint - translationSDF);
    return (grid::computeGradient_Robust(localSamplePoint, sdf));
}

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
    const EigenQuaternion<T>* A_centerRotation;

    const grid::Grid<T, T>* B_sdf;
    const EigenVector3<T>* B_linearVel;
    const EigenVector3<T>* B_angularVel;
    const EigenVector3<T>* B_centerTranslation;
    const EigenQuaternion<T>* B_centerRotation;
};

template <typename T> struct TriangleAtTimeInfo
{
    const EigenVector3<T> A_p0;
    const EigenVector3<T> A_p1;
    const EigenVector3<T> A_p2;
};

template <typename T>
EigenVector3<T> getTriangleVertexPosAt(const EigenMatrix3<T>& R,
                                       const EigenVector3<T>& centerTranslation,
                                       const EigenVector3<T>& linVel,
                                       const EigenVector3<T>& vert, T dt)
{
    //TODO USE EQUATION 14!, DISCARD THIS PELASE!
    EigenVector3<T> diff = vert - centerTranslation;

    return (centerTranslation + (linVel + R * diff) * dt);
}

template <typename T>
TriangleAtTimeInfo<T>
getTriangleAtTime_OLD(T t, const RigidBodyInfo<T>& initialState)
{
    //Current center position due to linear motion
    //Note we may have to set centertranslation to 0,0,0, because our p's may already be the rotated object pose!
    //EigenVector3<T> currentCenter = *(initialState.A_centerTranslation) + *(initialState.A_linearVel) * t;

    /*    EigenVector3<T> currentCenter = *(initialState.A_linearVel) * t;

    T angle = (*(initialState.A_angularVel)).norm() * t;

    EigenQuaternion<T> rotation;
    if (angle > 1e-10)
    {
        EigenVector3<T> axis = (*(initialState.A_linearVel)).normalized();
        rotation = EigenQuaternion<T>(Eigen::AngleAxis<T>(angle, axis));
    }
    else { rotation = EigenQuaternion<T>::Identity(); }*/

    /*    EigenVector3<T> currentCenter = *(initialState.A_centerTranslation);
    EigenVector3<T> radiusp0 = *(initialState.A_p0) - currentCenter;
    EigenVector3<T> velocityp0
        = (*(initialState.A_linearVel))
        + ((*(initialState.A_angularVel))).cross(radiusp0);
    EigenVector3<T> radiusp1 = *(initialState.A_p1) - currentCenter;
    EigenVector3<T> velocityp1
        = (*(initialState.A_linearVel))
        + ((*(initialState.A_angularVel))).cross(radiusp1);
    EigenVector3<T> radiusp2 = *(initialState.A_p2) - currentCenter;
    EigenVector3<T> velocityp2
        = (*(initialState.A_linearVel))
        + ((*(initialState.A_angularVel))).cross(radiusp2);
    TriangleAtTimeInfo<T> tInfo{.A_p0 = *(initialState.A_p0) + velocityp0 * t,
                                .A_p1 = *(initialState.A_p1) + velocityp1 * t,
                                .A_p2 = *(initialState.A_p2) + velocityp2 * t};*/

    EigenVector3<T> v_world = *(initialState.A_linearVel);
    EigenVector3<T> omega_world = *(initialState.A_angularVel);
    EigenVector3<T> C0 = *(initialState.A_centerTranslation);
    EigenVector3<T> p0 = *(initialState.A_p0);
    EigenVector3<T> p1 = *(initialState.A_p1);
    EigenVector3<T> p2 = *(initialState.A_p2);

    EigenMatrix3<T> R;

    if (true)
    {
        T const radian = omega_world.norm();
        Eigen::Matrix<T, 3, 1> axis = omega_world.normalized();
        R = Eigen::AngleAxis<T>(radian, axis).toRotationMatrix();
    }
    /*    EigenVector3<T> Cnew = C0 + v_world * t;
    
    //Apply exact rigid transform to each vertex
    EigenVector3<T> r0 = p0 - C0;
    EigenVector3<T> r1 = p1 - C0;
    EigenVector3<T> r2 = p2 - C0;

    EigenVector3<T> tmp0 = (Cnew + R * r0);
    EigenVector3<T> tmp1 = (Cnew + R * r1);
    EigenVector3<T> tmp2 = (Cnew + R * r2);*/
    EigenVector3<T> p0_t = getTriangleVertexPosAt(R, C0, v_world, p0, t);
    EigenVector3<T> p1_t = getTriangleVertexPosAt(R, C0, v_world, p0, t);
    EigenVector3<T> p2_t = getTriangleVertexPosAt(R, C0, v_world, p0, t);
    TriangleAtTimeInfo<T> out{.A_p0 = p0_t, .A_p1 = p1_t, .A_p2 = p2_t};

    return out;

    // Rotate and translate each vertex
    //I think we can do * instead of _transformVector!
    /*    TriangleAtTimeInfo<T> tInfo{
        .A_p0 = currentCenter + rotation._transformVector(*(initialState.A_p0)),
        .A_p1 = currentCenter + rotation._transformVector(*(initialState.A_p1)),
        .A_p2
        = currentCenter + rotation._transformVector(*(initialState.A_p2))};*/
    /*    TriangleAtTimeInfo<T> tInfo{.A_p0 = currentCenter + *(initialState.A_p0),
                                .A_p1 = currentCenter + *(initialState.A_p1),
                                .A_p2 = currentCenter + *(initialState.A_p2)};*/
    //    return tInfo;
}

template <typename T>
EigenVector3<T> getTriangleVertexPosAt(const EigenVector3<T>& centerTranslation,
                                       const EigenVector3<T>& linVel,
                                       const EigenVector3<T>& angVel,
                                       const EigenVector3<T>& vert, T dt)
{
    //TODO USE EQUATION 14!, DISCARD THIS PELASE!
    EigenVector3<T> diff = vert - centerTranslation;

    return vert + ((linVel + angVel.cross(diff)) * dt);
}

template <typename T>
TriangleAtTimeInfo<T> getTriangleAtTime(T t,
                                        const RigidBodyInfo<T>& initialState)
{
    //Current center position due to linear motion
    //Note we may have to set centertranslation to 0,0,0, because our p's may already be the rotated object pose!
    //EigenVector3<T> currentCenter = *(initialState.A_centerTranslation) + *(initialState.A_linearVel) * t;

    /*    EigenVector3<T> currentCenter = *(initialState.A_linearVel) * t;

    T angle = (*(initialState.A_angularVel)).norm() * t;

    EigenQuaternion<T> rotation;
    if (angle > 1e-10)
    {
        EigenVector3<T> axis = (*(initialState.A_linearVel)).normalized();
        rotation = EigenQuaternion<T>(Eigen::AngleAxis<T>(angle, axis));
    }
    else { rotation = EigenQuaternion<T>::Identity(); }*/

    /*    EigenVector3<T> currentCenter = *(initialState.A_centerTranslation);
    EigenVector3<T> radiusp0 = *(initialState.A_p0) - currentCenter;
    EigenVector3<T> velocityp0
        = (*(initialState.A_linearVel))
        + ((*(initialState.A_angularVel))).cross(radiusp0);
    EigenVector3<T> radiusp1 = *(initialState.A_p1) - currentCenter;
    EigenVector3<T> velocityp1
        = (*(initialState.A_linearVel))
        + ((*(initialState.A_angularVel))).cross(radiusp1);
    EigenVector3<T> radiusp2 = *(initialState.A_p2) - currentCenter;
    EigenVector3<T> velocityp2
        = (*(initialState.A_linearVel))
        + ((*(initialState.A_angularVel))).cross(radiusp2);
    TriangleAtTimeInfo<T> tInfo{.A_p0 = *(initialState.A_p0) + velocityp0 * t,
                                .A_p1 = *(initialState.A_p1) + velocityp1 * t,
                                .A_p2 = *(initialState.A_p2) + velocityp2 * t};*/

    EigenVector3<T> v_world = *(initialState.A_linearVel);
    EigenVector3<T> omega_world = *(initialState.A_angularVel);
    EigenVector3<T> C0 = *(initialState.A_centerTranslation);
    EigenVector3<T> p0 = *(initialState.A_p0);
    EigenVector3<T> p1 = *(initialState.A_p1);
    EigenVector3<T> p2 = *(initialState.A_p2);

    /*    EigenVector3<T> Cnew = C0 + v_world * t;
    
    //Apply exact rigid transform to each vertex
    EigenVector3<T> r0 = p0 - C0;
    EigenVector3<T> r1 = p1 - C0;
    EigenVector3<T> r2 = p2 - C0;

    EigenVector3<T> tmp0 = (Cnew + R * r0);
    EigenVector3<T> tmp1 = (Cnew + R * r1);
    EigenVector3<T> tmp2 = (Cnew + R * r2);*/
    EigenVector3<T> p0_t
        = getTriangleVertexPosAt(C0, v_world, omega_world, p0, t);
    EigenVector3<T> p1_t
        = getTriangleVertexPosAt(C0, v_world, omega_world, p1, t);
    EigenVector3<T> p2_t
        = getTriangleVertexPosAt(C0, v_world, omega_world, p2, t);
    TriangleAtTimeInfo<T> out{.A_p0 = p0_t, .A_p1 = p1_t, .A_p2 = p2_t};

    return out;

    // Rotate and translate each vertex
    //I think we can do * instead of _transformVector!
    /*    TriangleAtTimeInfo<T> tInfo{
        .A_p0 = currentCenter + rotation._transformVector(*(initialState.A_p0)),
        .A_p1 = currentCenter + rotation._transformVector(*(initialState.A_p1)),
        .A_p2
        = currentCenter + rotation._transformVector(*(initialState.A_p2))};*/
    /*    TriangleAtTimeInfo<T> tInfo{.A_p0 = currentCenter + *(initialState.A_p0),
                                .A_p1 = currentCenter + *(initialState.A_p1),
                                .A_p2 = currentCenter + *(initialState.A_p2)};*/
    //    return tInfo;
}

template <typename T>
EigenVector3<T> getVelocityAtPoint(const RigidBodyInfo<T>& initialState,
                                   const EigenVector3<T>& point, T t)
{
    //Get the center g of rotation at time t
    /*    EigenVector3<T> currentCenter
        = *(initialState.A_centerTranslation) + *(initialState.A_linearVel) * t;

    //Paper's formula: v_t = v_g + ω_g × (x_t - g)
    EigenVector3<T> radius = point - currentCenter;
    EigenVector3<T> velocity = *(initialState.A_linearVel)
                             + (*(initialState.A_angularVel)).cross(radius);

    return velocity;*/
    /*EigenVector3<T> currentCenter = *(initialState.A_centerTranslation);
    //EigenVector3<T> currentCenter = (EigenVector3<T>(0, 0, 0)) + *(initialState.A_linearVel) * t;

    //Paper's formula: v_t = v_g + ω_g × (x_t - g)
    EigenVector3<T> radius = point - currentCenter;
    EigenVector3<T> velocity = (*(initialState.A_linearVel))
                             + ((*(initialState.A_angularVel))).cross(radius);
    
    

    return velocity * t;*/
    EigenVector3<T> v_world = *(initialState.A_linearVel);
    EigenVector3<T> omega_world = *(initialState.A_angularVel);
    EigenVector3<T> C0 = *(initialState.A_centerTranslation);
    EigenVector3<T> p0 = *(initialState.A_p0);
    EigenVector3<T> p1 = *(initialState.A_p1);
    EigenVector3<T> p2 = *(initialState.A_p2);
    EigenVector3<T> Cnew = v_world;

    /*    EigenMatrix3<T> R;

    if (true)
    {
        T const radian = omega_world.norm() * t;
        Eigen::Matrix<T, 3, 1> axis = omega_world.normalized();
        R = Eigen::AngleAxis<T>(radian, axis).toRotationMatrix();
    }

    EigenVector3<T> r0 = point - C0;

    EigenVector3<T> tmp0 = (Cnew + R * r0);
    return tmp0;*/
    return v_world + omega_world.cross(point - C0);
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
    return std::abs<T>(valueAtProjection(*(params.grid), x,
                                         *(info.B_centerTranslation),
                                         *(info.B_centerRotation)));
}

template <typename T>
T SignedDistanceAtTime(const RigidBodyInfo<T>& info,
                       const DistanceAtTimeParams<T>& params, T t)
{
    EigenVector3<T> x
        = BarycentricInterpolate(params.u, params.v, params.w, t, info);
    return (valueAtProjection(*(params.grid), x, *(info.B_centerTranslation),
                              *(info.B_centerRotation)));
}

template <typename T> struct DistanceAtPointParams
{
    const grid::Grid<T, T>* grid;
};

template <typename T>
T SignedDistanceAtPoint(const RigidBodyInfo<T>& info,
                        const DistanceAtPointParams<T>& params,
                        EigenVector3<T> x)
{
    return (valueAtProjection(*(params.grid), x, *(info.B_centerTranslation),
                              *(info.B_centerRotation)));
}

template <typename T, typename F>
T GSSMinimize_WHAT_MODIFIED(T lstart, T lend, F func,
                            const DistanceAtTimeParams<T>& params,
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
    uint16_t it = 0;
    while ((l3 - l0) > tol(l1 + l2) && it < 8)
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

template <typename T>
T SignedDistanceAtPoint_MODIFIED(const RigidBodyInfo<T>& info,
                                 const DistanceAtPointParams<T>& params,
                                 T llambda, EigenVector3<T> pstart,
                                 EigenVector3<T> pend)
{
    EigenVector3<T> point = pstart + (pend - pstart) * llambda;
    return (valueAtProjection(*(params.grid), point,
                              *(info.B_centerTranslation),
                              *(info.B_centerRotation)));
}

template <typename T, typename F>
EigenVector3<T> GSSMinimize_WHAT_MODIFIED(
    T lstart, T lend, EigenVector3<T> pstart, EigenVector3<T> pend, F func,
    const DistanceAtPointParams<T>& params, const RigidBodyInfo<T>& info)
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
    T f0 = func(info, params, l0, pstart, pend);
    T f1 = func(info, params, l1, pstart, pend);
    T f2 = func(info, params, l2, pstart, pend);
    T f3 = func(info, params, l3, pstart, pend);
    uint16_t it = 0;
    while ((l3 - l0) > tol(l1 + l2) && it < 8)
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
}

template <typename T, typename F>
T GSSMinimize_WHAT(T lstart, T lend, F func,
                   const DistanceAtTimeParams<T>& params,
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

template <typename T, typename F>
EigenVector3<T> GSSMinimize_WHAT(EigenVector3<T> lstart, EigenVector3<T> lend,
                                 F func, const DistanceAtPointParams<T>& params,
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

    EigenVector3<T> l0 = lstart;
    EigenVector3<T> l1 = lerp(lstart, lend, alpha1);
    EigenVector3<T> l2 = lerp(lstart, lend, alpha2);
    EigenVector3<T> l3 = lend;
    T f0 = func(info, params, l0);
    T f1 = func(info, params, l1);
    T f2 = func(info, params, l2);
    T f3 = func(info, params, l3);
    T tol = 1e-5;
    //TODO NOT IMPLEMENTED CORRECTLY
    //while ((l3 - l0) <= tol((l1 + l2)))
    uint16_t it = 0;
    while ((l3 - l0).norm() <= ((l1 + l2).norm() * tol) && it < 8)
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
    const T phi = (1 + std::sqrt(5)) / 2;
    const T invphi = 1 / phi;

    // initial interior points
    T c = b - (b - a) / phi;
    T d = a + (b - a) / phi;

    T f_c = func(info, params, c);
    T f_d = func(info, params, d);
    int it = 0;

    while ((b - a) > tolerance && it < 8)
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

template <typename T>
EigenVector3<T> lerp(EigenVector3<T>& x0, EigenVector3<T>& x1, T alpha)
{
    EigenVector3<T> val = EigenVector3<T>(std::lerp<T>(x0.x(), x1.x(), alpha),
                                          std::lerp<T>(x0.y(), x1.y(), alpha),
                                          std::lerp<T>(x0.z(), x1.z(), alpha));
    return val;
}

template <typename T, typename F>
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
}

template <typename T>
void computeBarycentricCoordinates(const EigenVector3<T>& p0,
                                   const EigenVector3<T>& p1,
                                   const EigenVector3<T>& p2,
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

template <typename T> T sign(T val)
{
    T res;
    if (val < 0) { res = T(-1); }
    else if (val > 0) { res = T(1); }
    else { res = val; }
    return res;
}

template <typename T> void projectToTriangle(T& u, T& v, T& w)
{
    //Ensure barycentric coordinates are valid
    u = std::max<T>(0.0, std::min<T>(1.0, u));
    v = std::max<T>(0.0, std::min<T>(1.0, v));
    w = std::max<T>(0.0, std::min<T>(1.0, w));

    //Normalize to sum to 1
    T sum = u + v + w;
    if (sum > 0)
    {
        u /= sum;
        v /= sum;
        w /= sum;
    }
    else { u = v = w = 1.0 / 3.0; }
}

template <typename T>
T FrankWolfeGSS(T tstart, T tend, const RigidBodyInfo<T>& initialState/*const EigenVector3<T>& p0,
                const EigenVector3<T>& p1, EigenVector3<T>& p2,
                const grid::Grid<T, T>& grid,*/)
{
    T t1 = tstart;
    T ti = t1;
    T tip1 = std::numeric_limits<T>::max();
    //For now if there is no velocity, it means our object is stationary. Thus it can never hit the other object. We can thus ignore it.
    if ((*(initialState.A_linearVel)).norm() < 0.0000001f) { return tend; }

    //TODO  Compute the barycentric coordinates 𝑢, 𝑣, 𝑤 of the starting iterate.
    //Very importantly: For now we incorrectly assume velocity only for the triangles, not the SDF!

    //Okay this is a very important detail, we have our initial vi is not equation
    // from the paper, it is actually the linear velocity of our triangle rigid body!
    // Thus we must see that vti and vi is different and vti uses getvelocityatpoint!
    EigenVector3<T> vi = *(initialState.A_linearVel);
    EigenVector3<T> gradP0 = gradientAtProjection(
        *(initialState.A_p0), *(initialState.B_sdf),
        *(initialState.B_centerTranslation), *(initialState.B_centerRotation));
    EigenVector3<T> gradP1 = gradientAtProjection(
        *(initialState.A_p1), *(initialState.B_sdf),
        *(initialState.B_centerTranslation), *(initialState.B_centerRotation));
    EigenVector3<T> gradP2 = gradientAtProjection(
        *(initialState.A_p2), *(initialState.B_sdf),
        *(initialState.B_centerTranslation), *(initialState.B_centerRotation));
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
                                              .grid = initialState.B_sdf,
                                              .p0 = initialState.A_p0,
                                              .p1 = initialState.A_p1,
                                              .p2 = initialState.A_p2};

    float eps = 1e-6;
    size_t maxIterations = 32u;
    DistanceAtPointParams distanceAtPointParams{.grid = initialState.B_sdf};
    for (size_t i = 0; i < maxIterations; ++i)
    {
        EigenVector3<T> xti = BarycentricInterpolate(u, v, w, ti, initialState);
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
            tend = std::min<T>(ti, tend);
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
                    ti, tend, SignedDistanceAtTime<T>, distanceAtTimeParams,
                    initialState);
            }
        }
        // Solve spatial sub-problem
        //Note, xtip1 = x_{t_{i+1}}
        EigenVector3<T> xtip1
            = BarycentricInterpolate(u, v, w, tip1, initialState);
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
            *(initialState.A_angularVel), *(initialState.A_p0), tip1);
        EigenVector3<T> p1_at_ti = getTriangleVertexPosAt(
            *(initialState.A_centerTranslation), *(initialState.A_linearVel),
            *(initialState.A_angularVel), *(initialState.A_p1), tip1);
        EigenVector3<T> p2_at_ti = getTriangleVertexPosAt(
            *(initialState.A_centerTranslation), *(initialState.A_linearVel),
            *(initialState.A_angularVel), *(initialState.A_p2), tip1);
        /*        T p0Min = (*(initialState.A_p0)).dot(gradPhixtip1);
        T p1Min = (*(initialState.A_p1)).dot(gradPhixtip1);
        T p2Min = (*(initialState.A_p2)).dot(gradPhixtip1);

        if (phixtip1 <= 0) { tend = std::min<T>(tip1, tend); }
        EigenVector3<T> si;
        //Computing support vertex pi
        if (p0Min <= p1Min && p0Min <= p2Min) { si = *(initialState.A_p0); }
        else if (p1Min <= p2Min && p1Min <= p0Min)
        {
            si = *(initialState.A_p1);
        }
        else { si = *(initialState.A_p2); }*/
        T p0Min = (p0_at_ti).dot(gradPhixtip1);
        T p1Min = (p1_at_ti).dot(gradPhixtip1);
        T p2Min = (p2_at_ti).dot(gradPhixtip1);

        if (phixtip1 <= 0) { tend = std::min<T>(tip1, tend); }
        EigenVector3<T> si;
        //Computing support vertex pi
        if (p0Min <= p1Min && p0Min <= p2Min) { si = p0_at_ti; }
        else if (p1Min <= p2Min && p1Min <= p0Min) { si = p1_at_ti; }
        else { si = p2_at_ti; }
        xtip1 = GSSMinimize_WHAT_MODIFIED(T(0), T(1), xtip1, si,
                                          SignedDistanceAtPoint_MODIFIED<T>,
                                          distanceAtPointParams, initialState);
        //TODO: Update barycentric coordinates 𝑢, 𝑣, 𝑤 using x®𝑡𝑖+1
        computeBarycentricCoordinates(*(initialState.A_p0),
                                      *(initialState.A_p1),
                                      *(initialState.A_p2), xtip1, u, v, w);
        projectToTriangle(u, v, w);

        if (std::abs(tip1 - ti) <= eps
            && (std::abs(xtip1.x() - xti.x()) <= eps
                && std::abs(xtip1.y() - xti.y()) <= eps
                && std::abs(xtip1.z() - xti.z()) <= eps)
            && (phixtip1 <= eps))
        {
            break;
        }
        ti = tip1;
        xti = xtip1;
    }
    std::cerr << "ENDED UP WITH tip1 = " << tip1 << " and ti = " << ti << "\n";
    if (ti <= 0.00000001)
    {
        //THis code forces debug breakpoint, DELETE LATER when I figure out why TOI=0
        T a = 0;
        T b = a + ti;
        std::cerr << b << ";";
    }
    return std::min<T>(tip1, ti);
}
} // namespace grid
#endif // GRID_CCD_GOLDEN_SECTION_SEARCH_HPP
