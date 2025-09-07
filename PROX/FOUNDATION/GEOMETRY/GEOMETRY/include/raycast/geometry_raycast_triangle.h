#ifndef GEOMETRY_COMPUTE_RAYCAST_TRIANGLE_H
#define GEOMETRY_COMPUTE_RAYCAST_TRIANGLE_H

#include <types/geometry_triangle.h>
#include <types/geometry_plane.h>
#include <types/geometry_ray.h>
#include <raycast/geometry_raycast_plane.h>

#include <cmath>
#include <vector>

namespace geometry
{

template <typename T>
inline bool compute_raycast_triangle(RayEigen<T> const& ray,
                                     Triangle<T> const& triangle,
                                     EigenVector3<T>& hit, T& length,
                                     bool const& only_front_face = false)
{

    length = std::numeric_limits<T>::max();
    hit = EigenVector3<T>(0, 0, 0);

    Plane<T> const& plane = make_plane(triangle);

    EigenVector3<T> p = EigenVector3<T>(0, 0, 0);
    T t     = std::numeric_limits<T>::max();

    bool const hit_plane = compute_raycast_plane(ray, plane, p, t, only_front_face );

    if (!hit_plane)
      return false;

    bool const inside = inside_triangle( p, triangle );

    if(!inside)
      return false;

    hit    = p;
    length = t;

    return true;
  }


}// namespace geometry

// GEOMETRY_COMPUTE_RAYCAST_TRIANGLE_H
#endif
