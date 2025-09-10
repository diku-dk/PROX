#ifndef PROCEDURAL_TYPES_H
#define PROCEDURAL_TYPES_H

#include <cstdlib>
#include <eigenhelperall.h>

namespace procedural
{
        /**
         * A Geometry Handle.
     * A geometry handle points to a body space geometry stored
     * inside a simulator (identified by its gid = geometry indentifier)
     * together with mass properties.
     *
     * To make it more convenient to accurately place geometry in a local
     * frame of reference the instance keep track of the body frame to
     * model frame transformation.
     *
     * Hence, when procedural generated scenes need to place geometry in the
     * world they simply specify the model frames placement with respect to a
     * given local reference frame.
     *
     * In many cases the local reference frame would be identical to the world
     * coordinate system. However, in some cases when building scenes in a
     * recursive manner it is convenient to move a whole collection of the
     * scene by simply changing the local frame to world coordinate frame
     * transformation.
         */
template <typename T> class GeometryHandleEigen
{

public:
    T m_m;    ///< Total mass
    T m_Ixx;  ///< Body frame inertia tensor
    T m_Iyy;  ///< Body frame inertia tensor
    T m_Izz;  ///< Body frame inertia tensor
    EigenVector3<T> m_T;    ///< Body to model transformation
    EigenQuaternion<T> m_Q;    ///< Body to model transformation
    size_t m_gid; ///< unique identifier of the geometry in the physics engine

public:
    GeometryHandleEigen()
        : m_m(0)
        , m_Ixx(0)
        , m_Iyy(0)
        , m_Izz(0)
        , m_T()
        , m_Q()
        , m_gid()
    {
    }

    GeometryHandleEigen(T const& m, T const& Ixx, T const& Iyy, T const& Izz,
                        const EigenVector3<T>& t, const EigenQuaternion<T>& q,
                        size_t const& g)
        : m_m(m)
        , m_Ixx(Ixx)
        , m_Iyy(Iyy)
        , m_Izz(Izz)
        , m_T(t)
        , m_Q(q)
        , m_gid(g)
    {
    }

public:
    /**
       * Get Translation
       *
       * @return   body to model frame translation
       */
    const EigenVector3<T>& Tb2m() const { return m_T; }

    /**
       * Get Rotation.
       *
       * @return   body to model frame rotation
       */
    const EigenQuaternion<T>& Qb2m() const { return m_Q; }
        };

        template <typename T> class MaterialInfo
        {
        public:
            size_t m_stone_mid;
            size_t m_ground_mid;
            size_t m_cannonball_mid;

            T m_stone_density;
            T m_cannonball_density;
            T m_ground_density;
        };

        } // end of namespace procedural

// PROCEDURAL_TYPES_H
#endif
