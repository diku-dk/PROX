#ifndef GEOMETRY_COMPUTE_RAYCAST_DOP_H
#define GEOMETRY_COMPUTE_RAYCAST_DOP_H

#include <types/geometry_ray.h>
#include <types/geometry_plane.h>
#include <types/geometry_dop.h>

#include <geometry_inside.h>

#include <eigenhelperall.h>

#include <cmath>
#include <cassert>

namespace geometry
{

template <typename T, size_t K>
inline bool compute_raycast_dop(RayEigen<T> const& ray, DOP<T, K> const& dop,
                                EigenVector3<T>& point, T& length,
                                const T& threshold)
{

    assert(threshold >= 0
           || !"compute_raycast_dop(): threhold must be non-negative");

    size_t const N = K / 2;
    DirectionTableEigen<T, N> const D = DirectionTableEigenHelper<T, N>::make();

    length = std::numeric_limits<T>::max();
    point = EigenVector3<T>(0, 0, 0);

    const EigenVector3<T>& p = ray.origin();
    const EigenVector3<T>& r = ray.direction();

    T t_min = std::numeric_limits<T>::lowest();
    T t_max = std::numeric_limits<T>::max();

    for (unsigned int k = 0u; k < N; ++k)
    {
        const EigenVector3<T> d0 = -D(k);
        T const w0 = -dop(k).lower();
        T const f0 = dot(d0, r);

        const EigenVector3<T> d1 = D(k);
        T const w1 = dop(k).upper();
        T const f1 = dot(d1, r);

        if (f0 == 0) continue;

        T t_near = (w0 - dot(d0, p)) / f0;
        T t_far = (w1 - dot(d1, p)) / f1;

        assert(is_number(t_near) || !"compute_raycast_dop() NaN");
        assert(is_finite(t_near) || !"compute_raycast_dop() Inf");
        assert(is_number(t_far) || !"compute_raycast_dop() NaN");
        assert(is_finite(t_far) || !"compute_raycast_dop() Inf");

        if (t_near > t_far)
        {
            T const tmp = t_far;

            t_far = t_near;
            t_near = tmp;
        }

        t_min = (t_near > t_min) ? t_near : t_min;
        t_max = (t_far < t_max) ? t_far : t_max;
    }

    if (t_min > 0 && t_min <= t_max)
    {
        const EigenVector3<T> q = p + r * t_min;

        if (outside_dop(q, dop, threshold)) return false;

        length = t_min;
        point = q;

        assert(is_number(length) || !"compute_raycast_dop() NaN");
        assert(is_finite(length) || !"compute_raycast_dop() Inf");

        assert(is_number(point(0)) || !"compute_raycast_dop() NaN");
        assert(is_finite(point(0)) || !"compute_raycast_dop() Inf");
        assert(is_number(point(1)) || !"compute_raycast_dop() NaN");
        assert(is_finite(point(1)) || !"compute_raycast_dop() Inf");
        assert(is_number(point(2)) || !"compute_raycast_dop() NaN");
        assert(is_finite(point(2)) || !"compute_raycast_dop() Inf");

        return true;
    }

    return false;
}

template <typename T, size_t K>
inline bool compute_raycast_dop(RayEigen<T> const& ray, DOP<T, K> const& dop,
                                EigenVector3<T>& point, T& length)
{
    return compute_raycast_dop<T, K>(ray, dop, point, length, 0);
}

} //namespace geometry

// GEOMETRY_COMPUTE_RAYCAST_DOP_H
#endif
