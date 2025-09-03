#ifndef GEOMETRY_CONVERT_H
#define GEOMETRY_CONVERT_H

#include <types/geometry_dop.h>
#include <types/geometry_sphere.h>
#include <types/geometry_aabb.h>
#include <types/geometry_cylinder.h>
#include <types/geometry_capsule.h>

namespace geometry
{

  template<size_t K,typename T>
  inline DOP<T,K> convert( Sphere<T> const & S )
  {
    size_t              const N = K/2;
    DirectionTableEigen<T,N> const D = DirectionTableEigenHelper<T,N>::make();

    DOP<T, K> kdop;

    for(size_t k =  0u; k < N; ++k)
    {
        T const o = dot( D(k), (S.center()) );

      kdop(k).lower() = o - S.radius();
      kdop(k).upper() = o + S.radius();
    }

    return kdop;
  }

  template<size_t K,typename T>
  inline DOP<T,K> convert( AABBEigen<T> const & aabb )
  {
    size_t              const N = K/2;
    DirectionTableEigen<T,N> const D = DirectionTableEigenHelper<T,N>::make();

    DOP<T, K> kdop;

    for(size_t k =  0u; k < N; ++k)
    {
      T const l = dot(  D(k), aabb.min() );
      T const u = dot(  D(k), aabb.max() );

      kdop(k).lower() = l;
      kdop(k).upper() = u;
    }

    return kdop;
  }

  template<typename T>
  inline CapsuleEigen<T> convert( CylinderEigen<T> const & cylinder )
  {

    T const radius = cylinder.radius();

    const EigenVector3<T> point0 = cylinder.center() - cylinder.axis()*cylinder.height()*0.5f;
    const EigenVector3<T> point1 = cylinder.center() + cylinder.axis()*cylinder.height()*0.5f;

    return make_capsule(radius, point0, point1);
  }


}// namespace geometry

// GEOMETRY_CONVERT_H
#endif
