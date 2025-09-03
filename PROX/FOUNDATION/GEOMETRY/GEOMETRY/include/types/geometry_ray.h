#ifndef GEOMETRY_RAY_H
#define GEOMETRY_RAY_H

#include "tiny_math_types.h"
namespace geometry
{

  template<typename V>
  class Ray
  {
  public:
      using T = typename V::real_type;
      using VT = typename V::value_traits;

  protected:
      V m_origin;
      V m_direction;

  public:

    V const & origin()    const { return this->m_origin; }
    V       & origin()          { return this->m_origin; }

    V const & direction() const { return this->m_direction; }

    void set_direction(V const & direction)
    {
      this->m_direction = unit(direction);
    }

  public:

    Ray()
    : m_origin( V::zero() )
    , m_direction(  V::k() )
    {}

    Ray(V const & origin, V const & direction)
    : m_origin(origin)
    , m_direction( unit( direction ) )
    {}

    Ray(Ray const & ray)
    : m_origin(ray.m_origin)
    , m_direction(ray.m_direction)
    {}

    Ray & operator=(Ray const & ray)
    {
      if(this!=&ray)
      {
        this->m_origin    = ray.m_origin;
        this->m_direction = ray.m_direction;
      }
      return *this;
    }

  };

  template<typename V>
  inline Ray<V> make_ray(V const & origin, V const & direction)
  {
    return Ray<V>(origin, direction);
  }

}// namespace geometry

  namespace geometry
  {

  template<typename T>
  class RayEigen
  {
  protected:
      EigenVector3<T> m_origin;
      EigenVector3<T> m_direction;

  public:

      const EigenVector3<T>& origin()    const { return this->m_origin; }
      EigenVector3<T>& origin()          { return this->m_origin; }

      const EigenVector3<T>& direction() const { return this->m_direction; }

      void set_direction(const EigenVector3<T>& direction)
      {
          this->m_direction = (direction).normalized();
      }

  public:

      RayEigen()
          : m_origin( EigenVector3<T>(0,0,0) )
          , m_direction(  EigenVector3<T>(0,0,1) )
      {}

      RayEigen(const EigenVector3<T>& origin, const EigenVector3<T>& direction)
          : m_origin(origin)
          , m_direction( ( direction ).normalized() )
      {}

      RayEigen(RayEigen const & ray)
          : m_origin(ray.m_origin)
          , m_direction(ray.m_direction)
      {}

      RayEigen & operator=(RayEigen const & ray)
      {
          if(this!=&ray)
          {
              this->m_origin    = ray.m_origin;
              this->m_direction = ray.m_direction;
          }
          return *this;
      }

  };

  template<typename T>
  inline RayEigen<T> make_ray(const EigenVector3<T>& origin, const EigenVector3<T>& direction)
  {
      return RayEigen<T>(origin, direction);
  }

  template<typename V>
  RayEigen<typename V::real_type> convertRayToEigen(const Ray<V> ray)
  {
      RayEigen<typename V::real_type> newRay;
      newRay.set_direction(toEigen(ray.direction()));
      newRay.origin() = toEigen(ray.origin());
      return newRay;
  }
  }// namespace geometry



// GEOMETRY_RAY_H
#endif
