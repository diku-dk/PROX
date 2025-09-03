#ifndef GEOMETRY_POINT_H
#define GEOMETRY_POINT_H

#include <types/geometry_support_mapping.h>

namespace geometry
{

  template<typename T>
  class Point
    : public geometry::SupportMapping<T>
  {
  protected:
      EigenVector3<T> m_coordinates;

  public:

    const EigenVector3<T>& coord() const   { return this->m_coordinates; }
    EigenVector3<T>& coord()         { return this->m_coordinates; }

  public:

    Point()
    : m_coordinates(EigenVector3<T>(0,0,0) )
    {}

    Point(const EigenVector3<T>&  coordinates )
    : m_coordinates( coordinates )
    {}

    Point(Point const & p)
    {
      *this = p;
    }

    Point & operator=(Point const & p)
    {
      if( this != &p)
      {
        this->m_coordinates = p.m_coordinates;
      }
      return *this;
    }

  public:

    EigenVector3<T> get_support_point(EigenVector3<T> dir) const override
    {
        return (m_coordinates);
    }


    T get_scale() const
    {
      return 0;
    }
  };

} //namespace geometry

// GEOMETRY_POINT_H
#endif
