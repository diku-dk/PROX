#ifndef GRID_CCD_GOLDEN_SECTION_SEARCH_HPP
#define GRID_CCD_GOLDEN_SECTION_SEARCH_HPP

#include <grid_grid.h>
#include <grid_local_optimization.hpp>
#include <chrono>

namespace grid
{
/*enum FuncType
{
    SignedDistanceAtTime,
    UnsignedDistanceAtTime,
    SignedDistanceAtPoint
};

template <typename T> T someFunction(T val) { return val + val; }*/
template <typename T>
inline void barycentric(const EigenVector3<T>& x1, const EigenVector3<T>& x2,
                        const EigenVector3<T>& x3, const EigenVector3<T>& p,
                        T& w1, T& w2, T& w3)
{
    using std::sqrt;

    // Before computing the barycentric coordinates, we 'normalize' the triangle to ensure it is not a sliver
    auto b1 = x2 - x1;
    auto b2 = x3 - x1;
    auto b3 = ((b1).normalized()).cross((b2).normalized());

    const EigenMatrix3<T> basis{
        {b1(0), b2(0), b3(0)},
        {b1(1), b2(1), b3(1)},
        {b1(2), b2(2), b3(2)}
    };

    /*      M const basis = M::make(  b1(0), b2(0), b3(0)
                              , b1(1), b2(1), b3(1)
                              , b1(2), b2(2), b3(2));*/

    // q is the point p transformed to the isoparametric system
    const EigenVector3<T> q = basis.inverse() * (p - x1);

    w1 = 1 - q[0] - q[1];
    w2 = q[0];
    w3 = q[1];

    assert(is_number(w1) || !"barycentric(): NaN encountered");
    assert(is_number(w2) || !"barycentric(): NaN encountered");
    assert(is_number(w3) || !"barycentric(): NaN encountered");
}

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
    /*    contactPoint
        = transformRotation * contactPoint + transformTranslation;*/
    return (rotationSDF * grid::computeGradient_Working(localSamplePoint, sdf));
}

template <typename T> T tol(T val)
{
    T tol = 1e-8;
    //TODO NOT IMPLEMENTED CORRECTLY
    return tol * val;
}

template <typename T> struct RigidBodyInfo
{
    //For rigid body A:
    EigenVector3<T> A_p0;
    EigenVector3<T> A_p1;
    EigenVector3<T> A_p2;
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
    EigenVector3<T> A_p0;
    EigenVector3<T> A_p1;
    EigenVector3<T> A_p2;
};

template <typename T>
EigenVector3<T> getTriangleVertexPosAt_2(
    const EigenMatrix3<T>& R, const EigenVector3<T>& centerTranslation,
    const EigenVector3<T>& linVel, const EigenVector3<T>& vert, T dt)
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
    EigenVector3<T> p0 = (initialState.A_p0);
    EigenVector3<T> p1 = (initialState.A_p1);
    EigenVector3<T> p2 = (initialState.A_p2);

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
    EigenVector3<T> p1_t = getTriangleVertexPosAt(R, C0, v_world, p1, t);
    EigenVector3<T> p2_t = getTriangleVertexPosAt(R, C0, v_world, p2, t);
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
    /*EigenVector3<T> diff = vert - (centerTranslation + linVel * dt);

    return vert + ((linVel + angVel.cross(diff)) * dt);*/
    EigenVector3<T> diff = vert - (centerTranslation);

    return vert + ((linVel + angVel.cross(diff)) * dt);
}

template <typename T>
EigenVector3<T> getTriangleVertexInSDFLocalAtTime(
    const EigenVector3<T>& triangleCenter,
    const EigenVector3<T>& triangleLinVel,
    const EigenVector3<T>& triangleAngVel, const EigenVector3<T>& sdfCenter,
    const EigenVector3<T>& sdfLinVel, const EigenVector3<T>& sdfAngVel,
    const EigenVector3<T>& vert, T t)
{
    // 1. Compute triangle vertex position in world space at time t

    // Vector from triangle center to vertex in initial configuration
    EigenVector3<T> diffFromTriangleCenter = vert - triangleCenter;

    // Apply rotation to this vector
    T triangleAngle = triangleAngVel.norm() * t;
    EigenMatrix3<T> R_triangle;
    if (triangleAngle > 1e-8)
    {
        EigenVector3<T> axis = triangleAngVel.normalized();
        R_triangle
            = Eigen::AngleAxis<T>(triangleAngle, axis).toRotationMatrix();
    }
    else { R_triangle = EigenMatrix3<T>::Identity(); }

    // Final triangle vertex position in world space
    EigenVector3<T> trianglePosAtT = triangleCenter + triangleLinVel * t
                                   + R_triangle * diffFromTriangleCenter;

    // 2. Compute SDF transformation at time t

    // SDF center position at time t
    EigenVector3<T> sdfCenterAtT = sdfCenter + sdfLinVel * t;

    // SDF rotation at time t
    T sdfAngle = sdfAngVel.norm() * t;
    EigenMatrix3<T> R_sdf;
    if (sdfAngle > 1e-8)
    {
        EigenVector3<T> axis = sdfAngVel.normalized();
        R_sdf = Eigen::AngleAxis<T>(sdfAngle, axis).toRotationMatrix();
    }
    else { R_sdf = EigenMatrix3<T>::Identity(); }

    // 3. Transform triangle position to SDF local space
    // This is: R_sdf^T * (world_point - sdf_center)
    return R_sdf.transpose() * (trianglePosAtT - sdfCenterAtT);
}

template <typename T>
EigenVector3<T> getTriangleVertexPosAt_HMMM(
    const EigenVector3<T>& centerTranslationA, const EigenVector3<T>& linVel,
    const EigenVector3<T>& angVel, const EigenVector3<T>& centerTranslationB,
    const EigenVector3<T>& SDFlinVel, const EigenVector3<T>& SDFangVel,
    const EigenVector3<T>& vert, T dt)
{
    return getTriangleVertexInSDFLocalAtTime(centerTranslationA, linVel, angVel,
                                             centerTranslationB, SDFlinVel,
                                             SDFangVel, vert, dt);
    //TODO USE EQUATION 14!, DISCARD THIS PELASE!
    /*EigenVector3<T> diff = vert - (centerTranslation + linVel * dt);

    return vert + ((linVel + angVel.cross(diff)) * dt);*/
    EigenVector3<T> diffA = vert - (centerTranslationA);

    //EigenVector3<T> vi = ((linVel + angVel.cross(diffA)) * dt);
    EigenVector3<T> tmpAngVelA = angVel;
    T angleA = tmpAngVelA.norm();
    EigenMatrix3<T> RA;
    if (angleA > 1e-8)
    {
        EigenVector3<T> axisA = tmpAngVelA.normalized();
        RA = Eigen::AngleAxis<T>(angleA, axisA).toRotationMatrix();
    }
    else { RA = EigenMatrix3<T>::Identity(); }
    EigenVector3<T> vi = ((linVel + (RA) * (diffA)) * dt);

    EigenVector3<T> diffB = vert - (centerTranslationB);

    //EigenVector3<T> SDFVi = ((SDFlinVel + SDFangVel.cross(diffB)) * dt);
    EigenVector3<T> tmpAngVelB = SDFangVel;
    T angleB = tmpAngVelB.norm();
    EigenMatrix3<T> RB;
    if (angleB > 1e-8)
    {
        EigenVector3<T> axisB = tmpAngVelB.normalized();
        RB = Eigen::AngleAxis<T>(angleB, axisB).toRotationMatrix();
    }
    else { RB = EigenMatrix3<T>::Identity(); }
    EigenVector3<T> SDFVi = ((SDFlinVel + (RB) * (diffB)) * dt);

    EigenVector3<T> v_relative = vi - SDFVi;
    return vert + (v_relative);
}

template <typename T>
EigenVector3<T> getTriangleVertexPosAt(
    const EigenVector3<T>& centerTranslationA, const EigenVector3<T>& linVel,
    const EigenVector3<T>& angVel, const EigenVector3<T>& centerTranslationB,
    const EigenVector3<T>& SDFlinVel, const EigenVector3<T>& SDFangVel,
    const EigenVector3<T>& vert, T dt)
{
    //TODO USE EQUATION 14!, DISCARD THIS PELASE!
    /*EigenVector3<T> diff = vert - (centerTranslation + linVel * dt);

    return vert + ((linVel + angVel.cross(diff)) * dt);*/
    EigenVector3<T> diffA = vert - (centerTranslationA);

    //EigenVector3<T> vi = ((linVel + angVel.cross(diffA)) * dt);

    EigenVector3<T> vi = ((linVel + (angVel).cross(diffA)) * dt);

    EigenVector3<T> diffB = vert - (centerTranslationB);

    //EigenVector3<T> SDFVi = ((SDFlinVel + SDFangVel.cross(diffB)) * dt);

    EigenVector3<T> SDFVi = ((SDFlinVel + (SDFangVel.cross((diffB)))) * dt);

    EigenVector3<T> v_relative = vi - SDFVi;
    return vert + (v_relative);
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
    EigenVector3<T> CSDF = *(initialState.B_centerTranslation);
    EigenVector3<T> vSDF = *(initialState.B_linearVel);
    EigenVector3<T> omegaSDF = *(initialState.B_angularVel);
    EigenVector3<T> p0 = (initialState.A_p0);
    EigenVector3<T> p1 = (initialState.A_p1);
    EigenVector3<T> p2 = (initialState.A_p2);

    /*    EigenVector3<T> Cnew = C0 + v_world * t;
    
    //Apply exact rigid transform to each vertex
    EigenVector3<T> r0 = p0 - C0;
    EigenVector3<T> r1 = p1 - C0;
    EigenVector3<T> r2 = p2 - C0;

    EigenVector3<T> tmp0 = (Cnew + R * r0);
    EigenVector3<T> tmp1 = (Cnew + R * r1);
    EigenVector3<T> tmp2 = (Cnew + R * r2);*/
    EigenVector3<T> p0_t = getTriangleVertexPosAt(C0, v_world, omega_world,
                                                  CSDF, vSDF, omegaSDF, p0, t);
    EigenVector3<T> p1_t = getTriangleVertexPosAt(C0, v_world, omega_world,
                                                  CSDF, vSDF, omegaSDF, p1, t);
    EigenVector3<T> p2_t = getTriangleVertexPosAt(C0, v_world, omega_world,
                                                  CSDF, vSDF, omegaSDF, p2, t);
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
EigenVector3<T> getVelocityAtPoint2(const RigidBodyInfo<T>& initialState,
                                    const EigenVector3<T>& point, T t)
{
    EigenVector3<T> v_world = *(initialState.A_linearVel);
    EigenVector3<T> omega_world = *(initialState.A_angularVel);
    EigenVector3<T> C0 = *(initialState.A_centerTranslation);

    EigenVector3<T> vSDF = *(initialState.B_linearVel);
    EigenVector3<T> omegaSDF = *(initialState.B_angularVel);
    EigenVector3<T> CSDF = *(initialState.B_centerTranslation);

    // For triangle: use Equation 14 directly
    EigenVector3<T> diffA = point - C0;
    EigenVector3<T> v_triangle = v_world + omega_world.cross(diffA);

    // For SDF: also use Equation 14
    EigenVector3<T> diffB = point - CSDF;
    EigenVector3<T> v_sdf = vSDF + omegaSDF.cross(diffB);

    // Relative velocity
    return v_triangle - v_sdf;
}

template <typename T>
EigenVector3<T> getVelocityAtPoint_HMMM(const RigidBodyInfo<T>& initialState,
                                        const EigenVector3<T>& point, T t)
{
    //getVelocityAtPoint_GREAT(initialState, point, t);
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
    EigenVector3<T> CSDF = *(initialState.B_centerTranslation);
    EigenVector3<T> vSDF = *(initialState.B_linearVel);
    EigenVector3<T> omegaSDF = *(initialState.B_angularVel);
    EigenVector3<T> p0 = (initialState.A_p0);
    EigenVector3<T> p1 = (initialState.A_p1);
    EigenVector3<T> p2 = (initialState.A_p2);
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
    /*    EigenVector3<T> vti = (v_world + omega_world.cross(point - C0));
    EigenVector3<T> vtiSDF = (vSDF + omegaSDF.cross(point - CSDF));*/
    T angleA = omega_world.norm();
    EigenMatrix3<T> RA;
    if (angleA > 1e-8)
    {
        EigenVector3<T> axisA = omega_world.normalized();
        RA = Eigen::AngleAxis<T>(angleA, axisA).toRotationMatrix();
    }
    else { RA = EigenMatrix3<T>::Identity(); }

    T angleB = omegaSDF.norm();
    EigenMatrix3<T> RB;
    if (angleB > 1e-8)
    {
        EigenVector3<T> axisB = omegaSDF.normalized();
        RB = Eigen::AngleAxis<T>(angleB, axisB).toRotationMatrix();
    }
    else { RB = EigenMatrix3<T>::Identity(); }
    EigenVector3<T> vti = (v_world + (RA) * (point - C0));
    EigenVector3<T> vtiSDF = (vSDF + (RB) * (point - CSDF));

    return vti - vtiSDF;

    //return (v_world + omega_world.cross(point - C0));
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
    EigenVector3<T> CSDF = *(initialState.B_centerTranslation);
    EigenVector3<T> vSDF = *(initialState.B_linearVel);
    EigenVector3<T> omegaSDF = *(initialState.B_angularVel);
    EigenVector3<T> p0 = (initialState.A_p0);
    EigenVector3<T> p1 = (initialState.A_p1);
    EigenVector3<T> p2 = (initialState.A_p2);
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
    /*    EigenVector3<T> vti = (v_world + omega_world.cross(point - C0));
    EigenVector3<T> vtiSDF = (vSDF + omegaSDF.cross(point - CSDF));*/

    EigenVector3<T> vti = (v_world + (omega_world).cross(point - C0));
    EigenVector3<T> vtiSDF = (vSDF + (omegaSDF).cross(point - CSDF));

    return vti - vtiSDF;

    //return (v_world + omega_world.cross(point - C0));
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
    T u;
    T v;
    T w;
    const grid::Grid<T, T>* grid;
    const EigenVector3<T> p0;
    const EigenVector3<T> p1;
    const EigenVector3<T> p2;
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

template <typename T>
EigenVector3<T> lerp(EigenVector3<T>& x0, EigenVector3<T>& x1, T alpha)
{
    EigenVector3<T> val = EigenVector3<T>(std::lerp<T>(x0.x(), x1.x(), alpha),
                                          std::lerp<T>(x0.y(), x1.y(), alpha),
                                          std::lerp<T>(x0.z(), x1.z(), alpha));
    return val;
}

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
    //return FrankWolfeGSSBisection(tstart, tend, initialState);

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

template <typename T, typename F>
T bisectionRootFind(T a, T b, F func, const DistanceAtTimeParams<T>& params,
                    const RigidBodyInfo<T>& info)

{
    T fa = func(info, params, a);
    T fb = func(info, params, b);
    T tol = 1e-8;
    int max_iter = 50;

    //Check if we already have a root at boundaries
    if (std::abs(fa) < tol) return a;
    if (std::abs(fb) < tol) return b;

    //Ensure the interval brackets a root
    if (fa * fb > 0)
    {
        //No sign change - cannot use bisection
        //Return midpoint as fallback
        return (a + b) / 2;
    }

    T c = a;
    for (int i = 0; i < max_iter; i++)
    {
        //Find midpoint
        c = (a + b) / 2;
        T fc = func(info, params, c);

        //Check for convergence
        if (std::abs(fc) < tol || (b - a) / 2 < tol) { return c; }

        //Update interval
        if (fa * fc < 0)
        {
            b = c;
            fb = fc;
        }
        else
        {
            a = c;
            fa = fc;
        }
    }

    return c;
}

template <typename T>
T FrankWolfeGSSBisection(T tstart, T tend, const RigidBodyInfo<T>& initialState/*const EigenVector3<T>& p0,
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
            /*tip1 = GSSMinimize_WHAT_MODIFIED(
                tstart, ti, UnsignedDistanceAtTime<T>, distanceAtTimeParams,
                initialState);*/
            tip1 = bisectionRootFind(tstart, ti, UnsignedDistanceAtTime<T>,
                                     distanceAtTimeParams, initialState);
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
                /*tip1 = GSSMinimize_WHAT_MODIFIED(
                    tstart, ti, SignedDistanceAtTime<T>, distanceAtTimeParams,
                    initialState);*/
                tip1 = bisectionRootFind(tstart, ti, SignedDistanceAtTime<T>,
                                         distanceAtTimeParams, initialState);
            }
            else
            {
                /*tip1 = GSSMinimize_WHAT_MODIFIED(
                    ti, end, SignedDistanceAtTime<T>, distanceAtTimeParams,
                    initialState);*/
                tip1 = bisectionRootFind(ti, end, SignedDistanceAtTime<T>,
                                         distanceAtTimeParams, initialState);
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
T FrankWolfeGSSSimple(T tstart, T tend, const RigidBodyInfo<T>& initialState)
{

    T t0 = tstart;
    T ti = t0;
    T end = tend;
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

    T tip1;
    T eps = 1e-6;
    size_t maxIterations = 32u;
    for (size_t i = 0; i < maxIterations; ++i)
    {
        EigenVector3<T> xti
            = BarycentricInterpolate<T>(u, v, w, ti, initialState);
        T phixti = valueAtProjection<T>(*(initialState.B_sdf), xti,
                                        *(initialState.B_centerTranslation),
                                        *(initialState.B_centerRotation));
        EigenVector3<T> gradPhixti = gradientAtProjection(
            xti, *(initialState.B_sdf), *(initialState.B_centerTranslation),
            *(initialState.B_centerRotation));
        if (phixti <= 0) { end = std::min<T>(end, ti); }

        EigenVector3<T> p0_at_ti = getTriangleVertexPosAt(
            *(initialState.A_centerTranslation), *(initialState.A_linearVel),
            *(initialState.A_angularVel), p0s, ti);
        EigenVector3<T> p1_at_ti = getTriangleVertexPosAt(
            *(initialState.A_centerTranslation), *(initialState.A_linearVel),
            *(initialState.A_angularVel), p1s, ti);
        EigenVector3<T> p2_at_ti = getTriangleVertexPosAt(
            *(initialState.A_centerTranslation), *(initialState.A_linearVel),
            *(initialState.A_angularVel), p2s, ti);
        EigenVector3<T> si;
        //Computing support vertex pi

        T p0Mins = (p0_at_ti).dot(gradPhixti);
        T p1Mins = (p1_at_ti).dot(gradPhixti);
        T p2Mins = (p2_at_ti).dot(gradPhixti);

        if (p0Mins < p1Mins && p0Mins < p2Mins) { si = p0_at_ti; }
        else if (p1Mins < p2Mins && p1Mins < p0Mins) { si = p1_at_ti; }
        else { si = p2_at_ti; }
        T di;
        if (phixti <= 0) { di = T(-1); }
        else
        {
            di = -sign<T>(
                gradPhixti.dot(getVelocityAtPoint(initialState, xti, ti)));
        }
        T alpha = T(2) / (T(i + 2));
        T dt = alpha * (end - t0) * di;
        tip1 = std::max<T>(t0, std::min<T>(end, ti + dt));
        //        T xtip1 = proj()
        EigenVector3<T> p0_at_tip1 = getTriangleVertexPosAt(
            *(initialState.A_centerTranslation), *(initialState.A_linearVel),
            *(initialState.A_angularVel), p0s, tip1);
        EigenVector3<T> p1_at_tip1 = getTriangleVertexPosAt(
            *(initialState.A_centerTranslation), *(initialState.A_linearVel),
            *(initialState.A_angularVel), p1s, tip1);
        EigenVector3<T> p2_at_tip1 = getTriangleVertexPosAt(
            *(initialState.A_centerTranslation), *(initialState.A_linearVel),
            *(initialState.A_angularVel), p2s, tip1);
        EigenVector3<T> xtip1
            = BarycentricInterpolate(u, v, w, tip1, initialState);
        computeBarycentricCoordinates(p0_at_tip1, p1_at_tip1, p2_at_tip1, xtip1,
                                      u, v, w);
        if (std::abs(tip1 - ti) <= eps
            && (std::abs(xtip1.x() - xti.x()) <= eps
                && std::abs(xtip1.y() - xti.y()) <= eps
                && std::abs(xtip1.z() - xti.z()) <= eps))
        {
            break;
        }
        ti = tip1;
        xti = xtip1;
    }
    std::cerr << "ENDED UP WITH tip1 = " << tip1 << " and ti = " << ti << "\n";
    return ti;
}

template <typename T, typename F>
T TENARY_SEARCH(T lstart, T lend, F func, const DistanceAtTimeParams<T>& params,
                const RigidBodyInfo<T>& info)
{
    // Ternary search on [l0, l3]
    T l0 = lstart;
    T l3 = lend;

    // Tolerance and iteration cap
    const T tol = T(1e-8);
    const uint16_t maxIt = 200;

    // Evaluate endpoints once
    T f0 = func(info, params, l0);
    T f3 = func(info, params, l3);

    T l1, l2, f1, f2;
    uint16_t it = 0;

    while ((l3 - l0) > tol && it < maxIt)
    {
        // Two interior points at 1/3 and 2/3 of the interval
        l1 = l0 + (l3 - l0) / T(3);
        l2 = l0 + (l3 - l0) * T(2) / T(3);

        f1 = func(info, params, l1);
        f2 = func(info, params, l2);

        // If f1 < f2, minimum lies in [l0, l2]; otherwise in [l1, l3]
        if (f1 < f2)
        {
            l3 = l2;
            f3 = f2;
            // keep f0 unchanged
        }
        else
        {
            l0 = l1;
            f0 = f1;
            // keep f3 unchanged
        }
        ++it;
    }

    // Choose best among l0, mid, l3
    T lmid = (l0 + l3) / T(2);
    T fmid = func(info, params, lmid);

    // reuse computed f0 and f3 where available
    T lmin = l0;
    T fmin = f0;
    if (fmid < fmin)
    {
        lmin = lmid;
        fmin = fmid;
    }
    if (f3 < fmin)
    {
        lmin = l3;
        fmin = f3;
    }

    return lmin;
}

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
                const RigidBodyInfo<T>& info, T tol = 1e-8, int max_iter = 128)
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
T derivativeFunc(const RigidBodyInfo<T>& initialState,
                 const DistanceAtTimeParams<T>& params, T t)
{
    EigenVector3<T> xt
        = BarycentricInterpolate(params.u, params.v, params.w, t, initialState);
    EigenVector3<T> gradXt = computeGradient_Working(xt, *(params.grid));
    EigenVector3<T> vt = getVelocityAtPoint(initialState, xt, t) * t;
    T vel = gradXt.dot(vt);
    return vel;
};

template <typename T>
T FrankWolfeBRENT_BENCHMARK_TIME(T tstart, T tend,
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

    EigenVector3<T> p0s = ((initialState.A_p0)).eval();
    EigenVector3<T> p1s = ((initialState.A_p1)).eval();
    EigenVector3<T> p2s = ((initialState.A_p2)).eval();

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

    float eps = 1e-7;
    size_t maxIterations = 32u;
    //The higher the below iterations, the higher likelihood we will get no false negatives
    size_t hardStopMaxIterations = 1000u;
    size_t its = 0;
    DistanceAtPointParams distanceAtPointParams{.grid = initialState.B_sdf};
    EigenVector3<T> xtip1 = EigenVector3<T>(0, 0, 0);
    EigenVector3<T> xti = EigenVector3<T>(0, 0, 0);

    for (size_t i = 0; i < hardStopMaxIterations; ++i)
    {
        its += 1;
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
        ti = tip1;
        xti = xtip1;
    }

    xtiPoint = xtip1;
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

template <typename T>
T FrankWolfeGSSBisection_BENCHMARK_TIME(T tstart, T tend,
                                        const RigidBodyInfo<T>& initialState,
                                        EigenVector3<T>& xtiPoint,
                                        std::vector<T>& minimizerSteps)
{
    minimizerSteps.clear();
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
            /*tip1 = GSSMinimize_WHAT_MODIFIED(
                tstart, ti, UnsignedDistanceAtTime<T>, distanceAtTimeParams,
                initialState);*/
            auto gss_start = std::chrono::high_resolution_clock::now();
            tip1 = bisectionRootFind(tstart, ti, UnsignedDistanceAtTime<T>,
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
                tip1 = bisectionRootFind(tstart, ti, SignedDistanceAtTime<T>,
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
                tip1 = bisectionRootFind(ti, end, SignedDistanceAtTime<T>,
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
        //Computing support vertex pi
        if (p0Mins <= p1Mins && p0Mins <= p2Mins) { si = p0_at_ti; }
        else if (p1Mins <= p2Mins && p1Mins <= p0Mins) { si = p1_at_ti; }
        else { si = p2_at_ti; }

        auto gss_start = std::chrono::high_resolution_clock::now();
        xtip1 = GSSMinimize_WHAT(xtip1, si, SignedDistanceAtPoint<T>,
                                 distanceAtPointParams, initialState);
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

        //projectToTriangle(u, v, w);
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
#endif // GRID_CCD_GOLDEN_SECTION_SEARCH_HPP
