#ifndef GEOMETRY_LINE_H
#define GEOMETRY_LINE_H

#include <eigenhelperall.h>

namespace geometry
{

template<typename T>
  class Line
  {
  protected:

    EigenVector3<T> m_point;
    EigenVector3<T> m_direction;

  public:

    EigenVector3<T>& point()           { return this->m_point;      }
    EigenVector3<T>& direction()       { return this->m_direction;  }
    const EigenVector3<T>& point()     const { return this->m_point;      }
    const EigenVector3<T>& direction() const { return this->m_direction;  }

  public:

    Line()
    : m_point()
    , m_direction()
    {}

    Line(const EigenVector3<T>& point, const EigenVector3<T>& direction)
    : m_point(point)
        , m_direction( ( direction).normalized() )
    {}

    Line(Line const & line)
    {
      *this = line;
    }

    Line & operator=(Line const & line)
    {
      if(this != &line)
      {
        this->m_point      = line.m_point;
        this->m_direction  = line.m_direction;
      }
      return *this;
    }

  };

  struct FROM_POINTS {};
  struct FROM_DIRECTION {};

  template<typename T>
  inline Line<T> make_line(const EigenVector3<T>& a, const EigenVector3<T>& b, FROM_POINTS const & )
  {
    return Line<T>( a, b - a  );
  }

  template<typename T>
  inline Line<T> make_line(const EigenVector3<T>& point, const EigenVector3<T>& direction, FROM_DIRECTION const & )
  {
    return Line<T>( point, direction );
  }

  template<typename T>
  inline Line<T> make_line(const EigenVector3<T>& a, const EigenVector3<T>& b)
  {
    return make_line( a, b, FROM_POINTS() );
  }

}// namespace geometry

// GEOMETRY_LINE_H
#endif
