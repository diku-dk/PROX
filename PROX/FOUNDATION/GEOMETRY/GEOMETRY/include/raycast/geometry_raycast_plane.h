#ifndef GEOMETRY_COMPUTE_RAYCAST_PLANE_H
#define GEOMETRY_COMPUTE_RAYCAST_PLANE_H

#include <cmath>
#include <tiny.h>
#include <tiny_math_types.h>
#include <types/geometry_ray.h>

namespace geometry
{

/**
   *
   * @param only_front_face       Default value is false. If set to true then
   *                              a hit is only found if ray is hitting from the
   *                              front side of the plane. If hitting from
   *                              the back side then no hit is reported.
   */
template <typename T>
inline bool compute_raycast_plane(RayEigen<T> const& ray, Plane<T> const& plane, EigenVector3<T>& hit, T& length,
                                  bool const& only_front_face = false)
{
    length = std::numeric_limits<T>::max();
    hit = EigenVector3<T>(0, 0, 0);

    const EigenVector3<T> p = ray.origin();
    const EigenVector3<T> r = ray.direction();
    const EigenVector3<T> n = plane.normal();
    T const w = plane.offset();

    // Test if ray is parallel with plane
    if (fabs(dot(n, r)) <= std::numeric_limits<T>::epsilon() * 10) return false;

    if (only_front_face && dot(n, r) > 0) return false;

    T const test = (w - dot(n, p)) / dot(n, r);

    assert(is_number(test) || !"compute_ray_plane(): nan");
    assert(is_finite(test) || !"compute_ray_plane(): inf");

    // Test if the ray hit will be "behind" the ray origin (shooting backwards)
    if (test < 0) return false;

    length = test;
    hit = p + length * r;

    assert(is_number(length) || !"compute_ray_plane(): nan");
    assert(is_finite(length) || !"compute_ray_plane(): inf");
    assert(length >= 0 || !"compute_ray_plane(): non-positive");

    assert(is_number(hit(0)) || !"compute_ray_plane(): nan");
    assert(is_finite(hit(0)) || !"compute_ray_plane(): inf");
    assert(is_number(hit(1)) || !"compute_ray_plane(): nan");
    assert(is_finite(hit(1)) || !"compute_ray_plane(): inf");
    assert(is_number(hit(2)) || !"compute_ray_plane(): nan");
    assert(is_finite(hit(2)) || !"compute_ray_plane(): inf");

    return true;
}

} //namespace geometry

// GEOMETRY_COMPUTE_RAYCAST_PLANE_H
#endif
