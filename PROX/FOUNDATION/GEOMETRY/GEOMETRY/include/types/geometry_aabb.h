#ifndef GEOMETRY_AABB_H
#define GEOMETRY_AABB_H

#include "tiny_math_types.h"
namespace geometry
{

  template<typename V>
  class AABB
  {
  public:
      using T = typename V::real_type;
      using VT = typename V::value_traits;

  protected:
      V m_min_coord;
      V m_max_coord;

  public:

    V const & min() const { return this->m_min_coord; }
    V       & min()       { return this->m_min_coord; }

    V const & max() const { return this->m_max_coord; }
    V       & max()       { return this->m_max_coord; }

  public:

    AABB()
    : m_min_coord( V::make(-0.5f,-0.5f,-0.5f) )
    , m_max_coord( V::make( 0.5f, 0.5f, 0.5f) )
    {}

    AABB(V const & min_coord, V const & max_coord)
    : m_min_coord( min_coord )
    , m_max_coord( max_coord )
    {}

    AABB(AABB const & aabb)
    : m_min_coord(aabb.m_min_coord)
    , m_max_coord(aabb.m_max_coord)
    {}

    AABB & operator=(AABB const & aabb)
    {
      if(this!=&aabb)
      {
        this->m_min_coord = aabb.m_min_coord;
        this->m_max_coord = aabb.m_max_coord;
      }
      return *this;
    }

  };

  template<typename V>
  inline AABB<V> make_aabb(V const & min_coord, V const & max_coord)
  {
    return AABB<V>(min_coord, max_coord);
  }

}// namespace geometry

  namespace geometry
  {

  template<typename T>
  class AABBEigen
  {
  protected:
      EigenVector3<T> m_min_coord;
      EigenVector3<T> m_max_coord;

  public:

      const EigenVector3<T>& min() const { return this->m_min_coord; }
      EigenVector3<T>& min()       { return this->m_min_coord; }

      const EigenVector3<T>& max() const { return this->m_max_coord; }
      EigenVector3<T>& max()       { return this->m_max_coord; }

  public:

      AABBEigen()
          : m_min_coord( EigenVector3<T>(-0.5f,-0.5f,-0.5f) )
          , m_max_coord( EigenVector3<T>( 0.5f, 0.5f, 0.5f) )
      {}

      AABBEigen(const EigenVector3<T>& min_coord, const EigenVector3<T>& max_coord)
          : m_min_coord( min_coord )
          , m_max_coord( max_coord )
      {}

      AABBEigen(AABBEigen const & aabb)
          : m_min_coord(aabb.m_min_coord)
          , m_max_coord(aabb.m_max_coord)
      {}

      AABBEigen& operator=(AABBEigen const & aabb)
      {
          if(this!=&aabb)
          {
              this->m_min_coord = aabb.m_min_coord;
              this->m_max_coord = aabb.m_max_coord;
          }
          return *this;
      }

  };

  template<typename T>
  inline AABBEigen<T> make_aabb(const EigenVector3<T>& min_coord, const EigenVector3<T>& max_coord)
  {
      return AABBEigen<T>(min_coord, max_coord);
  }
  template<typename V>
  AABBEigen<typename V::real_type> convertAABBToEigen(const AABB<V> ray)
  {
      AABBEigen<typename V::real_type> newRay;
      newRay.m_min_coord() = (toEigen(ray.m_min_coord()));
      newRay.m_max_coord() = toEigen(ray.m_max_coord());
      return newRay;
  }

  }

// GEOMETRY_AABB_H
#endif
