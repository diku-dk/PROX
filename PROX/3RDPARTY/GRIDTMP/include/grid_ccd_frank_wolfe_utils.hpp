#ifndef GRID_CCD_FRANK_WOLFE_UTILS_HPP
#define GRID_CCD_FRANK_WOLFE_UTILS_HPP

#include <grid_grid.h>
#include <grid_local_optimization.hpp>
#include <chrono>

namespace grid
{

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

template <typename T>
EigenVector3<T> lerp(EigenVector3<T>& x0, EigenVector3<T>& x1, T alpha)
{
    EigenVector3<T> val = EigenVector3<T>(std::lerp<T>(x0.x(), x1.x(), alpha),
                                          std::lerp<T>(x0.y(), x1.y(), alpha),
                                          std::lerp<T>(x0.z(), x1.z(), alpha));
    return val;
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

} // namespace grid

#endif
