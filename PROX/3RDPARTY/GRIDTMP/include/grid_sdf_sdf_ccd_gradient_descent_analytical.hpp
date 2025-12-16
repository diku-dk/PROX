#ifndef GRID_SDF_SDF_CCD_GRADIENT_DESCENT_ANALYTICAL_HPP
#define GRID_SDF_SDF_CCD_GRADIENT_DESCENT_ANALYTICAL_HPP

#include "grid_ccd_frank_wolfe_utils.hpp"
#include <grid_grid.h>
#include <grid_sdf_sdf_voxelize.hpp>
#include <grid_sdf_sdf_ccd_gradient_descent.hpp>
#include <chrono>

namespace SDFSDFContactAnalytical
{
template <typename T> using Vec3 = Eigen::Matrix<T, 3, 1>;

template <typename T> constexpr T sphere_radius()
{
    return static_cast<T>(0.5);
}

// Signed distance to the sphere isosurface (positive outside, negative inside)
template <typename T> inline T sphere_signed_distance(const Vec3<T>& p)
{
    static_assert(std::is_floating_point<T>::value,
                  "T must be a floating point type");
    return p.norm() - sphere_radius<T>();
}

// Absolute distance to the isosurface (unsigned distance)
template <typename T> inline T sphere_distance_to_isosurface(const Vec3<T>& p)
{
    return std::abs(sphere_signed_distance<T>(p));
}

// Gradient of the signed distance (unit outward normal).
// At the origin (||p|| == 0) the gradient is undefined; we return a zero vector there.
template <typename T> inline Vec3<T> sphere_sdf_gradient(const Vec3<T>& p)
{
    static_assert(std::is_floating_point<T>::value,
                  "T must be a floating point type");
    const T len = p.norm();
    const T eps = std::numeric_limits<T>::epsilon()
                * static_cast<T>(10); // small tolerance
    if (len > eps)
    {
        return p / len; // normalized
    }
    else
    {
        return Vec3<T>::Zero(); // undefined at center; choose zero for safety
    }
}

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
EigenVector3<T> sphereSDFNormal(EigenVector3<T> x, EigenVector3<T> c)
{
    EigenVector3<T> v = x - c;
    T n = v.norm();
    if (n == 0.0) return EigenVector3<T>(1.0, 0.0, 0.0);
    return v / n;
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

template <typename T>
T valueAtProjectionAna(const grid::Grid<T, T>& sdf,
                       const EigenVector3<T>& samplePoint,
                       const EigenVector3<T>& translationSDF,
                       const EigenQuaternion<T>& rotationSDF)
{
    EigenVector3<T> localSamplePoint
        = rotationSDF.inverse() * (samplePoint - translationSDF);
    return sphere_signed_distance(localSamplePoint);
}

template <typename T>
EigenVector3<T> gradientAtProjectionAna(const EigenVector3<T>& samplePoint,
                                        const grid::Grid<T, T>& sdf,
                                        const EigenVector3<T>& translationSDF,
                                        const EigenQuaternion<T>& rotationSDF)
{
    EigenVector3<T> localSamplePoint
        = rotationSDF.inverse() * (samplePoint - translationSDF);
    /*    contactPoint
        = transformRotation * contactPoint + transformTranslation;*/
    return (rotationSDF * sphere_sdf_gradient(localSamplePoint));
}

template <typename T>
EigenVector3<T>
projectToSDFSurfaceLocalA(EigenVector3<T> pos,
                          const SDFSDFContact::SingleRigidBodyInfo<T>& sdfInfo)
{
    return pos
         - valueAtProjectionAna(*(sdfInfo.sdf), pos,
                                *(sdfInfo.A_centerTranslation),
                                *(sdfInfo.A_centerRotation))
               * gradientAtProjectionAna(pos, *(sdfInfo.sdf),
                                         *(sdfInfo.A_centerTranslation),
                                         *(sdfInfo.A_centerRotation))
                     .normalized();
}

//Apply the standard SDF transformation and then the time transformation
template <typename T>
EigenVector3<T>
gradientAtProjectionForBB(const EigenVector3<T>& samplePoint,
                          const grid::Grid<T, T>& sdf,
                          const EigenVector3<T>& translationSDF,
                          const EigenQuaternion<T>& rotationSDF,
                          const SDFSDFContact::currentSDFPose<T> pose)
{
    EigenVector3<T> localSamplePoint
        = pose.rotation.inverse()
        * ((rotationSDF.inverse() * (samplePoint - translationSDF))
           - pose.translation);
    /*    contactPoint
        = transformRotation * contactPoint + transformTranslation;*/
    return (pose.rotation * rotationSDF)
         * (sphere_sdf_gradient(localSamplePoint));
}

//Apply the standard SDF transformation and then the time transformation
template <typename T>
T valueAtProjectionForBB(const EigenVector3<T>& samplePoint,
                         const grid::Grid<T, T>& sdf,
                         const EigenVector3<T>& translationSDF,
                         const EigenQuaternion<T>& rotationSDF,
                         const SDFSDFContact::currentSDFPose<T> pose)
{
    EigenVector3<T> localSamplePoint
        = pose.rotation.inverse()
        * ((rotationSDF.inverse() * (samplePoint - translationSDF))
           - pose.translation);
    return (sphere_signed_distance(localSamplePoint));
}

template <typename T>
T SignedDistance(const SDFSDFContact::currentSDFPose<T>& poseA_dead,
                 const SDFSDFContact::currentSDFPose<T>& poseB_dead,
                 const SDFSDFContact::SingleRigidBodyInfo<T>& SDFA,
                 const SDFSDFContact::SingleRigidBodyInfo<T>& SDFB, T ti,
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
    SDFSDFContact::currentSDFPose<T> poseB = {translationB, rotationB};

    /*    return std::abs<T>(grid::valueAtProjection(
        SDFB, newPos, *(SDFB.A_centerTranslation), *(SDFB.A_centerRotation)));*/
    return std::abs<T>(valueAtProjectionForBB(newPos, *(SDFB.sdf),
                                              *(SDFB.A_centerTranslation),
                                              *(SDFB.A_centerRotation), poseB));
}

template <typename T>
EigenVector3<T> compute_effective_point_velocity_for_A_at(
    const EigenVector3<T>& world_x, const EigenVector3<T>& cA,
    const EigenVector3<T>& v_cA, const EigenVector3<T>& omegaA)
{
    // cross = (x - cA) cross omega
    EigenVector3<T> r = world_x - cA;
    EigenVector3<T> cross = r.cross(omegaA);

    // effective velocity for partial_t (x held fixed)
    // note the minus on v_cA
    EigenVector3<T> v_eff = cross - v_cA;
    return v_eff;
}

template <typename T>
EigenVector3<T> getVelocityAtPoint(const EigenVector3<T>& translationA,
                                   const EigenVector3<T>& xt,
                                   const EigenVector3<T>& omega,
                                   const EigenVector3<T>& linearVelocity)
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

    //v®𝑡 = v®𝑔 + 𝜔®𝑔 × (®x𝑡 − g®)
    return linearVelocity + omega.cross(xt - translationA);

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

    //return (v_world + omega_world.cross(point - C0));
}

template <typename T>
T getSDFSDFTOISingleVoxelCCD(const EigenVector3<T>& position,
                             const SDFSDFContact::SingleRigidBodyInfo<T>& SDFA,
                             const SDFSDFContact::SingleRigidBodyInfo<T>& SDFB,
                             T tstart, T tend, T bestDTFound,
                             bool& contactFound,
                             EigenVector3<T>& outContactPoint)
{
    //TODO: If velocity is zero for this SDF combination, we can return because nothing can move in time. Maybe return -1 for nothing? Otherwise return 1 for end
    if ((*(SDFA.A_angularVel)).norm() <= 1e-12
        && (*(SDFA.A_linearVel)).norm() <= 1e-12)
    {
        return tend;
    }

    //First our pos is in local coordinates in the sdf with absolutely no transformations
    // not even those applied to the SDF. Thus we get the pos in world coordinates at time t=0:
    EigenVector3<T> pos
        = *(SDFA.A_centerRotation) * position + *(SDFA.A_centerTranslation);
    T ti = tstart;
    T tip1 = -1.0;
    EigenVector3<T> xip1;
    EigenVector3<T> x_ti;
    T scale = 1000;
    uint32_t maxIterations = 100000 / uint32_t(scale);
    T stepSizeAlpha = 0.0001 * scale;
    T stepSizeAlphaT = 0.0005 * scale;
    T stepSizeAlphaOriginal = stepSizeAlpha;
    T eps = 1e-8;
    //x_ti = position;
    pos = projectToSDFSurfaceLocalA(pos, SDFA);
    bool penetration = false;
    uint32_t its;
    for (uint32_t i = 0; i < maxIterations; ++i)
    {
        x_ti = getVertexPosAtMat(*(SDFA.A_centerTranslation),
                                 *(SDFA.A_linearVel), *(SDFA.A_angularVel), pos,
                                 ti);
        //Assume pos is in world coordinates but not transformed. So it is in world coordinates compared to
        // SDF_A. So basically pos is always the point in world coordinates at time t=0.

        //Let x always be in local coordiantes. Then we need no transformation for x
        // in SDF A. x should only be transformed to B when we do gradient computation with B
        EigenVector3<T> translationB;
        EigenQuaternion<T> rotationB;
        EigenVector3<T> translationA;
        EigenQuaternion<T> rotationA;
        //Early exit -- no need to search more if we have a significantly earlier candidate!
        if (ti /* * 0.5*/ > bestDTFound || ti >= tend) { return ti; }

        //We do the same as before, but ONLY get the rotation and translation needed
        // to transform SDF B into time ti, such that we can make accurate quries for
        // x_ti.
        getTransformForBody(*(SDFB.A_centerTranslation), *(SDFB.A_linearVel),
                            *(SDFB.A_angularVel), ti, translationB, rotationB);
        getTransformForBody(*(SDFA.A_centerTranslation), *(SDFA.A_linearVel),
                            *(SDFA.A_angularVel), ti, translationA, rotationA);
        SDFSDFContact::currentSDFPose<T> poseB = {translationB, rotationB};
        SDFSDFContact::currentSDFPose<T> poseA = {translationA, rotationA};

        EigenVector3<T> gradA = gradientAtProjectionForBB(
            x_ti, *(SDFA.sdf), *(SDFA.A_centerTranslation),
            *(SDFA.A_centerRotation), poseA);
        EigenVector3<T> normA = gradA.normalized();

        EigenVector3<T> gradB = gradientAtProjectionForBB(
            x_ti, *(SDFB.sdf), *(SDFB.A_centerTranslation),
            *(SDFB.A_centerRotation), poseB);

        EigenVector3<T> vtiA
            = getVelocityAtPoint(*(SDFA.A_centerTranslation), x_ti,
                                 *(SDFA.A_angularVel), *(SDFA.A_linearVel));
        EigenVector3<T> vtiB
            = getVelocityAtPoint(*(SDFB.A_centerTranslation), x_ti,
                                 *(SDFB.A_angularVel), *(SDFB.A_linearVel));
        T At = gradA.dot(vtiA);
        T Bt = gradB.dot(vtiB);
        EigenVector4<T> gB
            = EigenVector4<T>(gradB.x(), gradB.y(), gradB.z(), Bt);
        EigenVector4<T> nA
            = EigenVector4<T>(gradA.x(), gradA.y(), gradA.z(), At);

        //EigenVector3<T> gradientDir = (gradB - (gradB.dot(normA)) * normA);
        T g_dot_n = gB.dot(nA);
        T norm_n2 = dot(nA, nA);
        EigenVector4<T> gradientDir = (gB - (g_dot_n / norm_n2) * nA);
        /*        EigenVector4<T> gradientDirNew = (gB-nA);
        gradientDir.x() = gradientDirNew.x();
        gradientDir.y() = gradientDirNew.y();
        gradientDir.z() = gradientDirNew.z();*/
        //gradientDir.w() = -gradientDirNew.w();
        //gradientDir = -gB;
        EigenVector3<T> dx = EigenVector3<T>(gradientDir.x(), gradientDir.y(),
                                             gradientDir.z());
        T oldPointPenetration = valueAtProjectionForBB(
            x_ti, *(SDFB.sdf), *(SDFB.A_centerTranslation),
            *(SDFB.A_centerRotation), poseB);
        T dt
            = gradientDir.w()
            * oldPointPenetration; //std::max<T>(oldPointPenetration, stepSizeAlphaT * 0.01);
        //T dt = gradientDir.w() * std::min<T>((oldPointPenetration), 0.1);
        //Below can be either negative or positive, it actually doesnt matter much?
        EigenVector4<T> p = -gradientDir;
        //stepSizeAlpha = backtracking_line_search<T>(tstart, tend, ti, x_ti, p, gradientDir, SDFA, SDFB, poseB);
        xip1 = x_ti - stepSizeAlpha * dx;

        /*T alphaStepSizeT = backtracking_line_search_time_only(
            tstart, tend, ti, x_ti, p, gradientDir, SDFA, SDFB, poseB);
        tip1 = std::clamp<T>(ti + alphaStepSizeT * dt, tstart, tend);*/
        tip1 = std::clamp<T>(ti + stepSizeAlphaT * dt, tstart, tend);
        //Check for convergence...
        //Maybe check tip?
        getTransformForBody(*(SDFB.A_centerTranslation), *(SDFB.A_linearVel),
                            *(SDFB.A_angularVel), tip1, translationB,
                            rotationB);
        poseB = {translationB, rotationB};

        EigenVector3<T> cools = reverseVertexPosAtMat(
            *(SDFA.A_centerTranslation), *(SDFA.A_linearVel),
            *(SDFA.A_angularVel), xip1, ti);
        EigenVector3<T> atA = projectToSDFSurfaceLocalA(cools, SDFA);

        cools = getVertexPosAtMat(*(SDFA.A_centerTranslation),
                                  *(SDFA.A_linearVel), *(SDFA.A_angularVel),
                                  atA, tip1);

        T newPointPenetration = valueAtProjectionForBB(
            cools, *(SDFB.sdf), *(SDFB.A_centerTranslation),
            *(SDFB.A_centerRotation), poseB);

        /*if (newPointPenetration <= -eps)
        {
            stepSizeAlpha *= 0.5;
            tip1 = ti;
        }

        else */
        //If absolutely no movement, break!
        /*if (std::abs<T>((x_ti - cools).norm()) <= eps
            && std::abs<T>(ti - tip1) <= eps)
        {
            break;
        }*/
        if (newPointPenetration <= 1e-7)
        {
            T f0 = oldPointPenetration;
            T f1 = newPointPenetration;
            T t0 = ti;
            T t1 = tip1;
            EigenVector3<T> p0 = x_ti;
            EigenVector3<T> p1 = cools;

            T den = f1 - f0;
            T res = std::numeric_limits<T>::max();
            if (std::abs<T>(den) < 1e-12)
            {
                //DO nothing
            }
            else { res = ((eps - f0) / den); }

            if (res != std::numeric_limits<T>::max() && T(0.0) <= res <= T(1.0))
            {
                EigenVector3<T> p_star = p0 + res * (p1 - p0);
                T realTOI = t0 + res * (t1 - t0);
                tip1 = realTOI;
                xip1 = p_star;
            }
            else { xip1 = cools; }
            /*xip1 = reverseVertexPosAtMat(*(SDFA.A_centerTranslation),
                                         *(SDFA.A_linearVel),
                                         *(SDFA.A_angularVel), xip1, ti);
            xip1 = projectToSDFSurfaceLocal(xip1, SDFA);
            xip1 = getVertexPosAtMat(*(SDFA.A_centerTranslation),
                                     *(SDFA.A_linearVel), *(SDFA.A_angularVel),
                                     xip1, tip1);*/
            penetration = true;
            break;
        }
        /*else
        {
            //RESET
            stepSizeAlpha = stepSizeAlphaOriginal;
        }*/

        //Now traverse back to SDF start pose, such that xtip now lies in the
        // SDFs pose at t=0!
        //Pretty sure we should use ti!
        /*xip1 = reverseVertexPosAtMat(*(SDFA.A_centerTranslation),
                                     *(SDFA.A_linearVel), *(SDFA.A_angularVel),
                                     xip1, ti);
        //After transform, remember to project back to local coordinates!
        xip1 = projectToSDFSurfaceLocal(xip1, SDFA);*/
        xip1 = atA;
        //Now set our new search start point to pos!
        pos = xip1;
        ti = tip1;
        its = i;
    }
    //std::cerr << "ITS: " << its << "\n";
    // std::cerr << tip1 << "\n";
    //If no penetration was ever found, we simply do not have a TOI.
    if (!penetration) tip1 = tend;
    contactFound = penetration;
    outContactPoint = xip1;

    return tip1;
}

template <typename T>
T getSDFSDFTOISingleVoxelCCD_BENCHMARK_Z(
    const EigenVector3<T>& position,
    const SDFSDFContact::SingleRigidBodyInfo<T>& SDFA,
    const SDFSDFContact::SingleRigidBodyInfo<T>& SDFB, T tstart, T tend,
    T bestDTFound, bool& contactFound, size_t& iterationsUsed,
    EigenVector3<T>& outContactPoint)
{
    //TODO: If velocity is zero for this SDF combination, we can return because nothing can move in time. Maybe return -1 for nothing? Otherwise return 1 for end
    if ((*(SDFA.A_angularVel)).norm() <= 1e-12
        && (*(SDFA.A_linearVel)).norm() <= 1e-12)
    {
        return tend;
    }

    //First our pos is in local coordinates in the sdf with absolutely no transformations
    // not even those applied to the SDF. Thus we get the pos in world coordinates at time t=0:
    EigenVector3<T> pos
        = *(SDFA.A_centerRotation) * position + *(SDFA.A_centerTranslation);
    T ti = tstart;
    T tip1 = -1.0;
    EigenVector3<T> xip1;
    EigenVector3<T> x_ti;
    T scale = 1000;
    uint32_t maxIterations = 100000 / uint32_t(scale);
    T stepSizeAlpha = 0.0001 * scale;
    T stepSizeAlphaT = 0.0005 * scale;
    T stepSizeAlphaOriginal = stepSizeAlpha;
    T eps = 1e-8;
    //x_ti = position;
    pos = projectToSDFSurfaceLocalA(pos, SDFA);
    bool penetration = false;
    iterationsUsed = 0;
    for (uint32_t i = 0; i < maxIterations; ++i)
    {
        x_ti = getVertexPosAtMat(*(SDFA.A_centerTranslation),
                                 *(SDFA.A_linearVel), *(SDFA.A_angularVel), pos,
                                 ti);
        //Assume pos is in world coordinates but not transformed. So it is in world coordinates compared to
        // SDF_A. So basically pos is always the point in world coordinates at time t=0.

        //Let x always be in local coordiantes. Then we need no transformation for x
        // in SDF A. x should only be transformed to B when we do gradient computation with B
        EigenVector3<T> translationB;
        EigenQuaternion<T> rotationB;
        EigenVector3<T> translationA;
        EigenQuaternion<T> rotationA;
        //Early exit -- no need to search more if we have a significantly earlier candidate!
        if (ti /* * 0.5*/ > bestDTFound || ti >= tend) { return ti; }

        //We do the same as before, but ONLY get the rotation and translation needed
        // to transform SDF B into time ti, such that we can make accurate quries for
        // x_ti.
        getTransformForBody(*(SDFB.A_centerTranslation), *(SDFB.A_linearVel),
                            *(SDFB.A_angularVel), ti, translationB, rotationB);
        getTransformForBody(*(SDFA.A_centerTranslation), *(SDFA.A_linearVel),
                            *(SDFA.A_angularVel), ti, translationA, rotationA);
        SDFSDFContact::currentSDFPose<T> poseB = {translationB, rotationB};
        SDFSDFContact::currentSDFPose<T> poseA = {translationA, rotationA};

        EigenVector3<T> gradA = gradientAtProjectionForBB<T>(
            x_ti, *(SDFA.sdf), *(SDFA.A_centerTranslation),
            *(SDFA.A_centerRotation), poseA);
        EigenVector3<T> normA = gradA.normalized();

        EigenVector3<T> gradB = gradientAtProjectionForBB<T>(
            x_ti, *(SDFB.sdf), *(SDFB.A_centerTranslation),
            *(SDFB.A_centerRotation), poseB);

        EigenVector3<T> vtiA
            = getVelocityAtPoint(*(SDFA.A_centerTranslation), x_ti,
                                 *(SDFA.A_angularVel), *(SDFA.A_linearVel));
        EigenVector3<T> vtiB
            = getVelocityAtPoint(*(SDFB.A_centerTranslation), x_ti,
                                 *(SDFB.A_angularVel), *(SDFB.A_linearVel));
        T At = gradA.dot(vtiA);
        T Bt = gradB.dot(vtiB);
        EigenVector4<T> gB
            = EigenVector4<T>(gradB.x(), gradB.y(), gradB.z(), Bt);
        EigenVector4<T> nA
            = EigenVector4<T>(gradA.x(), gradA.y(), gradA.z(), At);

        //EigenVector3<T> gradientDir = (gradB - (gradB.dot(normA)) * normA);
        T g_dot_n = gB.dot(nA);
        T norm_n2 = dot(nA, nA);
        EigenVector4<T> gradientDir = (gB - (g_dot_n / norm_n2) * nA);
        /*        EigenVector4<T> gradientDirNew = (gB-nA);
        gradientDir.x() = gradientDirNew.x();
        gradientDir.y() = gradientDirNew.y();
        gradientDir.z() = gradientDirNew.z();*/
        //gradientDir.w() = -gradientDirNew.w();
        //gradientDir = -gB;
        EigenVector3<T> dx = EigenVector3<T>(gradientDir.x(), gradientDir.y(),
                                             gradientDir.z());
        T oldPointPenetration = valueAtProjectionForBB(
            x_ti, *(SDFB.sdf), *(SDFB.A_centerTranslation),
            *(SDFB.A_centerRotation), poseB);
        T dt
            = gradientDir.w()
            * oldPointPenetration; // std::max<T>(oldPointPenetration,             stepSizeAlphaT * 0.01);
        //T dt = gradientDir.w() * std::min<T>((oldPointPenetration), 0.1);
        //Below can be either negative or positive, it actually doesnt matter much?
        EigenVector4<T> p = -gradientDir;
        //stepSizeAlpha = backtracking_line_search<T>(tstart, tend, ti, x_ti, p, gradientDir, SDFA, SDFB, poseB);
        xip1 = x_ti - stepSizeAlpha * dx;

        /*T alphaStepSizeT = backtracking_line_search_time_only(
            tstart, tend, ti, x_ti, p, gradientDir, SDFA, SDFB, poseB);
        tip1 = std::clamp<T>(ti + alphaStepSizeT * dt, tstart, tend);*/
        tip1 = std::clamp<T>(ti + stepSizeAlphaT * dt, tstart, tend);
        //Check for convergence...
        //Maybe check tip?
        getTransformForBody(*(SDFB.A_centerTranslation), *(SDFB.A_linearVel),
                            *(SDFB.A_angularVel), tip1, translationB,
                            rotationB);
        poseB = {translationB, rotationB};

        EigenVector3<T> cools = reverseVertexPosAtMat(
            *(SDFA.A_centerTranslation), *(SDFA.A_linearVel),
            *(SDFA.A_angularVel), xip1, ti);
        EigenVector3<T> atA = projectToSDFSurfaceLocalA(cools, SDFA);

        cools = getVertexPosAtMat(*(SDFA.A_centerTranslation),
                                  *(SDFA.A_linearVel), *(SDFA.A_angularVel),
                                  atA, tip1);

        T newPointPenetration = valueAtProjectionForBB(
            cools, *(SDFB.sdf), *(SDFB.A_centerTranslation),
            *(SDFB.A_centerRotation), poseB);

        /*if (newPointPenetration <= -eps)
        {
            stepSizeAlpha *= 0.5;
            tip1 = ti;
        }

        else */
        //If absolutely no movement, break!
        /*if (std::abs<T>((x_ti - cools).norm()) <= eps
            && std::abs<T>(ti - tip1) <= eps)
        {
            break;
        }*/
        if (newPointPenetration <= 1e-7)
        {
            T f0 = oldPointPenetration;
            T f1 = newPointPenetration;
            T t0 = ti;
            T t1 = tip1;
            EigenVector3<T> p0 = x_ti;
            EigenVector3<T> p1 = cools;

            T den = f1 - f0;
            T res = std::numeric_limits<T>::max();
            if (std::abs<T>(den) < 1e-12)
            {
                //DO nothing
            }
            else { res = ((eps - f0) / den); }

            if (res != std::numeric_limits<T>::max() && T(0.0) <= res <= T(1.0))
            {
                EigenVector3<T> p_star = p0 + res * (p1 - p0);
                T realTOI = t0 + res * (t1 - t0);
                tip1 = realTOI;
                xip1 = p_star;
            }
            else { xip1 = cools; }
            /*xip1 = reverseVertexPosAtMat(*(SDFA.A_centerTranslation),
                                         *(SDFA.A_linearVel),
                                         *(SDFA.A_angularVel), xip1, ti);
            xip1 = projectToSDFSurfaceLocal(xip1, SDFA);
            xip1 = getVertexPosAtMat(*(SDFA.A_centerTranslation),
                                     *(SDFA.A_linearVel), *(SDFA.A_angularVel),
                                     xip1, tip1);*/
            penetration = true;
            iterationsUsed = (i + 1);
            break;
        }
        /*else
        {
            //RESET
            stepSizeAlpha = stepSizeAlphaOriginal;
        }*/

        //Now traverse back to SDF start pose, such that xtip now lies in the
        // SDFs pose at t=0!
        //Pretty sure we should use ti!
        /*xip1 = reverseVertexPosAtMat(*(SDFA.A_centerTranslation),
                                     *(SDFA.A_linearVel), *(SDFA.A_angularVel),
                                     xip1, ti);
        //After transform, remember to project back to local coordinates!
        xip1 = projectToSDFSurfaceLocal(xip1, SDFA);*/
        xip1 = atA;
        //Now set our new search start point to pos!
        pos = xip1;
        ti = tip1;
        iterationsUsed = (i + 1);
    }
    //std::cerr << "ITS: " << its << "\n";
    // std::cerr << tip1 << "\n";
    //If no penetration was ever found, we simply do not have a TOI.
    if (!penetration) tip1 = tend;
    contactFound = penetration;
    outContactPoint = xip1;

    return tip1;
}

template <typename T>
T getSDFSDFTOISingleVoxel_BENCHMARK_P(
    const SDFSDFContact::SDFVoxel<T>& SDFBVoxel,
    const SDFSDFContact::SingleRigidBodyInfo<T>& SDFA,
    const SDFSDFContact::SingleRigidBodyInfo<T>& SDFB, T tstart, T tend,
    T& avgIterationsForVoxel, EigenVector3<T>& outContactPoints)
{
    T minTOI = std::numeric_limits<T>::max();
    size_t totalIterationsUsed = 0;
    for (size_t i = 0; i < SDFBVoxel.selected.size(); ++i)
    {
        EigenVector3<T> selectedPoint = SDFBVoxel.selected[i].pos;
        EigenVector3<T> outContactPoint;
        bool contactFound = false;
        size_t iterationsUsed;
        T toi = getSDFSDFTOISingleVoxelCCD_BENCHMARK_Z(
            selectedPoint, SDFA, SDFB, tstart, tend, minTOI, contactFound,
            iterationsUsed, outContactPoint);
        totalIterationsUsed += iterationsUsed;
        /*if (std::abs<T>(toi - minTOI) <= 1e-9)
        {
            outContactPoints.push_back(outContactPoint);
        }
        else */
        if (toi <= minTOI && contactFound)
        {
            /*outContactPoints.clear();
            outContactPoints.push_back(outContactPoint);*/
            outContactPoints = outContactPoint;
            minTOI = toi;
        }
    }
    avgIterationsForVoxel
        = T(totalIterationsUsed) / T(SDFBVoxel.selected.size());
    return minTOI;
}

template <typename T>
T getSDFSDFTOI_BENCHMARK(
    const std::vector<SDFSDFContact::SDFVoxel<T>>& SDFAVoxel,
    const std::vector<SDFSDFContact::SDFVoxel<T>>& SDFBVoxel,
    const SDFSDFContact::SingleRigidBodyInfo<T>& SDFA,
    const SDFSDFContact::SingleRigidBodyInfo<T>& SDFB, T tstart, T tend,
    T& avgIterationsForVoxelA, T& avgIterationsForVoxelB,
    EigenVector3<T>& outContactPoints)
{
    avgIterationsForVoxelA = 0;
    T minTOI = std::numeric_limits<T>::max();
    //Compute all voxels of SDF A
    for (size_t i = 0; i < SDFAVoxel.size(); ++i)
    {
        EigenVector3<T> newOutContactPoints;
        //Importantly to get correct projection first loop is SDFA, SDFB, second should swap!
        T currAvgIterationsForVoxel = 0;
        T toi = getSDFSDFTOISingleVoxel_BENCHMARK_P(
            SDFAVoxel[i], SDFA, SDFB, tstart, tend, currAvgIterationsForVoxel,
            newOutContactPoints);
        avgIterationsForVoxelA += currAvgIterationsForVoxel;
        if (minTOI >= toi)
        {
            minTOI = toi;
            //outContactPoints = std::move(newOutContactPoints);
            outContactPoints = newOutContactPoints;
        };
    }
    avgIterationsForVoxelA = avgIterationsForVoxelA / T(SDFAVoxel.size());

    avgIterationsForVoxelB = 0;
    //Compute all voxels of SDF B
    for (size_t i = 0; i < SDFBVoxel.size(); ++i)
    {
        EigenVector3<T> newOutContactPoints;
        //Importantly to get correct projection first loop is SDFA, SDFB, second should swap!
        T currAvgIterationsForVoxel = 0;
        T toi = getSDFSDFTOISingleVoxel_BENCHMARK_P(
            SDFBVoxel[i], SDFB, SDFA, tstart, tend, currAvgIterationsForVoxel,
            newOutContactPoints);
        avgIterationsForVoxelB += currAvgIterationsForVoxel;
        if (minTOI >= toi)
        {
            minTOI = toi;
            outContactPoints = newOutContactPoints;
        }
    }
    avgIterationsForVoxelB = avgIterationsForVoxelB / T(SDFBVoxel.size());
    if (minTOI == std::numeric_limits<T>::max()) { return tend; }
    return minTOI;
}

} // namespace SDFSDFContact

#endif // GRID_SDF_SDF_CCD_GRADIENT_DESCENT_ANALYTICAL_HPP
