#ifndef GEOMETRY_SPHERE_H
#define GEOMETRY_SPHERE_H

#include <types/geometry_support_mapping.h>

#include <tiny_is_number.h>
#include <tiny_is_finite.h>
#include <tiny_math_types.h>

#include <cmath> // needed for std::sqrt
#include <cassert>

namespace geometry
{

  template<typename T>
  class Sphere
    : public geometry::SupportMapping<T>

  {

  protected:
      EigenVector3<T> m_center;
      T m_radius;

  public:

    const EigenVector3<T>& center() const { return this->m_center; }
    const T radius() const { return this->m_radius; }
    EigenVector3<T>& center()       { return this->m_center; }
    T& radius()       { return this->m_radius; }

  public:

    Sphere()
          : m_center( EigenVector3<T>({0,0,0}) )
      , m_radius( 1 )
    {
    }

    Sphere(const EigenVector3<T>& center, const T radius )
      : m_center( center )
      , m_radius( radius )
    {
    }

    Sphere(Sphere const & sph)
    {
      *this = sph;
    }

    Sphere & operator=(Sphere const & sph)
    {
      if( this != &sph)
      {
        this->m_center = sph.m_center;
        this->m_radius = sph.m_radius;
      }
      return *this;
    }


  public:
      EigenVector3<T> get_support_point(EigenVector3<T> dir) const override
      {
          auto dirLen = dir.norm();
          assert(is_number(dirLen) && is_finite(dirLen));
          if (dirLen > 0) { return (m_radius / dirLen) * dir + (m_center); }
          return (m_center) + EigenVector3<T>{m_radius, 0, 0};
    }

    T get_scale() const
    {
      assert( is_number(this->m_radius)   || !"get_scale(): NAN encountered");
      assert( is_finite(this->m_radius)   || !"get_scale(): INF encountered");
      assert( this->m_radius > 0 || !"get_scale(): radius was non-positive");

      return 2 *  this->m_radius ;
    }

  };

  template<typename T>
  inline Sphere<T> make_sphere(const EigenVector3<T>& center, const T radius)
  {
    return Sphere<T>(center,radius);
  }

}// namespace geometry

//GEOMETRY_SPHERE_H
#endif
