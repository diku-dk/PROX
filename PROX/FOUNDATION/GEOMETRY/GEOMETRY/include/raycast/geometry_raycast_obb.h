#ifndef GEOMETRY_COMPUTE_RAYCAST_OBB_H
#define GEOMETRY_COMPUTE_RAYCAST_OBB_H

#include "geometry_transform.h"
#include <types/geometry_obb.h>
#include <types/geometry_ray.h>

#include <tiny_is_number.h>
#include <tiny_is_finite.h>

#include <cmath>
#include <vector>

namespace geometry
{

  template<typename T>
  inline bool compute_raycast_obb(
                                   RayEigen<T> const & ray
                                   , OBBEigen<T> const & obb
                                   , EigenVector3<T>& q
                                   , T& length
                                   )
  {

    q      = EigenVector3<T>(0,0,0);
    length = std::numeric_limits<T>::max();


    const EigenVector3<T> p_obb  = transform_to_obb( ray.origin(), obb, TRANSFORM_POINT()     );
    const EigenVector3<T> r_obb  = transform_to_obb( ray.direction(), obb, TRANSFORM_VECTOR() );

    const EigenVector3<T> aabb_mincoord = -obb.half_extent();
    const EigenVector3<T> aabb_maxcoord =  obb.half_extent();

    EigenVector3<T> q_obb;
    T t_obb;

    bool const did_hit = compute_raycast_aabb(p_obb, r_obb, aabb_mincoord, aabb_maxcoord, q_obb, t_obb);

    if(did_hit)
    {
      length = t_obb;
      q      = transform_from_obb(q_obb, obb, TRANSFORM_POINT() );
      return true;
    }

    return false;
  }


} //namespace geometry

// GEOMETRY_COMPUTE_RAYCAST_OBB_H
#endif
