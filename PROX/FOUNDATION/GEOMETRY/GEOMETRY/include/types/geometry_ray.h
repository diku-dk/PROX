#ifndef GEOMETRY_RAY_H
#define GEOMETRY_RAY_H

#include <eigenhelperall.h>

namespace geometry
{

template <typename T> class RayEigen
{
protected:
    EigenVector3<T> m_origin;
    EigenVector3<T> m_direction;

public:
    const EigenVector3<T>& origin() const { return this->m_origin; }
    EigenVector3<T>& origin() { return this->m_origin; }

    const EigenVector3<T>& direction() const { return this->m_direction; }

    void set_direction(const EigenVector3<T>& direction)
    {
        this->m_direction = (direction).normalized();
    }

public:
    RayEigen()
        : m_origin(EigenVector3<T>(0, 0, 0))
        , m_direction(EigenVector3<T>(0, 0, 1))
    {
    }

    RayEigen(const EigenVector3<T>& origin, const EigenVector3<T>& direction)
        : m_origin(origin)
        , m_direction((direction).normalized())
    {
    }

    RayEigen(RayEigen const& ray)
        : m_origin(ray.m_origin)
        , m_direction(ray.m_direction)
    {
    }

    RayEigen& operator=(RayEigen const& ray)
    {
        if (this != &ray)
        {
            this->m_origin = ray.m_origin;
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

  }// namespace geometry



// GEOMETRY_RAY_H
#endif
