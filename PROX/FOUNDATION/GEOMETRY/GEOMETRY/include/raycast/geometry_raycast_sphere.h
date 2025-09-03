#ifndef GEOMETRY_COMPUTE_RAYCAST_SPHERE_H
#define GEOMETRY_COMPUTE_RAYCAST_SPHERE_H

#include "tiny_vector_functions.h"
#include <cmath> // needed for std::sqrt
#include <tiny_is_finite.h>
#include <tiny_is_number.h>
#include <types/geometry_ray.h>
#include <types/geometry_sphere.h>

namespace geometry
{

template <typename T>
inline bool compute_raycast_sphere(RayEigen<T> const& ray, Sphere<T> const& sphere, EigenVector3<T>& hit_point,
                                   T& length)
{
    using std::sqrt;

    length = std::numeric_limits<T>::max();
    hit_point = EigenVector3<T>(0, 0, 0);

    const EigenVector3<T>& o = (sphere.center());
    T const& r = sphere.radius();
    const EigenVector3<T>& p = ray.origin();
    const EigenVector3<T>& u = ray.direction();

    T const a = dot(u, u);
    EigenVector3<T> tmp = p - o;
    T const b = dot(tmp, u) * 2;
    T const c = dot(tmp, tmp) - r * r;

    assert(a > 0 || !"compute_raycast_sphere(): a-coefficent must be positive");

    // Test that ray origin is out-side sphere
    if (c <= 0) { return false; }

    const auto D = b * b - 4 * a * c;

    if (D < 0) return false;

    T const sqrt_D = sqrt(D);

    T const r1 = (-b - sqrt_D) / (2 * a);
    T const r2 = (-b + sqrt_D) / (2 * a);

    assert(is_finite(r1) || !"compute_raycast_sphere(): r1 is inf");
    assert(is_number(r1) || !"compute_raycast_sphere(): r1 is nan");
    assert(is_finite(r2) || !"compute_raycast_sphere(): r2 is inf");
    assert(is_number(r2) || !"compute_raycast_sphere(): r2 is nan");

    assert(r1 <= r2 || !"compute_raycast_sphere(): r1 can not be larger than r2");

    if (r1 >= 0)
    {
        length = r1;
        hit_point = p + length * u;
        return true;
    }

    if (r2 >= 0)
    {
        length = r1;
        hit_point = p + length * u;
        return true;
    }

    return false;
}

} //namespace geometry

// GEOMETRY_COMPUTE_RAYCAST_SPHERE_H
#endif
