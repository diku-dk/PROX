#ifndef HYPER_MODIFIERS_H
#define HYPER_MODIFIERS_H

#include <hyper_body.h>

#include <tiny.h>
#include <mesh_array.h>

#include <cassert>
#include <vector>

namespace hyper
{

  template<typename MT>
  inline void rotate(
                Body<MT> & body
                , typename MT::vector3_type const & axis
                , typename MT::real_type const & angle
                , typename MT::vector3_type const & center
                )
  {
    typedef typename MT::vector3_type    V;
    typedef typename MT::quaternion_type Q;

    Q const q = Q::Ru(angle, axis);

    for(size_t k = 0u; k < body.m_X.size(); ++k)
    {
      V const p = V::make( body.m_X(k), body.m_Y(k), body.m_Z(k) );

      V const r = rotate(q,p-center) + center;

      body.m_X(k) = r(0);
      body.m_Y(k) = r(1);
      body.m_Z(k) = r(2);
    }
    body.m_x = MT::convert(body.m_X, body.m_Y, body.m_Z);
  }

  template<typename MT>
  inline void rotate(
                     Body<MT>  & body
                     , typename MT::vector3_type const & axis
                     , typename MT::real_type const & angle
                     )
  {
    typedef typename MT::vector3_type    V;

    rotate(body, axis, angle, V::zero() );
  }

  template<typename MT>
  inline void translate(
                        Body<MT> & body
                        , typename MT::vector3_type const & d
                        )
  {
    for(size_t k = 0u; k < body.m_X.size(); ++k)
    {
      body.m_X(k) = body.m_X(k) + d(0);
      body.m_Y(k) = body.m_Y(k) + d(1);
      body.m_Z(k) = body.m_Z(k) + d(2);
    }

    body.m_x = MT::convert(body.m_X, body.m_Y, body.m_Z);
  }

  template<typename MT>
  inline void scale(
                        Body<MT> & body
                        , typename MT::vector3_type const & s
                        )
  {
    for(size_t k = 0u; k < body.m_X.size(); ++k)
    {
      body.m_X(k) = body.m_X(k) * s(0);
      body.m_Y(k) = body.m_Y(k) * s(1);
      body.m_Z(k) = body.m_Z(k) * s(2);
    }

    body.m_x = MT::convert(body.m_X, body.m_Y, body.m_Z);
  }

  template<typename MT>
  inline void twist(
                    Body<MT> & body
                    , typename MT::vector3_type const & axis
                    , typename MT::real_type const & angle
                    )
  {
    typedef typename MT::real_type       T;
    typedef typename MT::vector3_type    V;
    typedef typename MT::value_traits    VT;
    typedef typename MT::quaternion_type Q;

    size_t const N = body.m_X.size();

    V const n = unit(axis);

    V center  = V::make( body.m_X(0), body.m_Y(0), body.m_Z(0) );
    T w_start = inner_prod(n, center );
    T w_end   = inner_prod(n, center );

    for(size_t k = 1u; k < N; ++k)
    {
      V const p = V::make( body.m_X(k), body.m_Y(k), body.m_Z(k) );
      T const w = inner_prod(n, p);

      w_start = w < w_start ? w : w_start;
      w_end   = w > w_end   ? w : w_end;

      center = center + p;
    }
    center = center / VT::numeric_cast( N );

    for(size_t k = 0u; k < N; ++k)
    {
      V const p = V::make( body.m_X(k), body.m_Y(k), body.m_Z(k) );
      T const w = inner_prod(n, p);
      T const a = angle*(w-w_start) / (w_end-w_start);
      Q const q = Q::Ru(a, n);
      V const o = center + w*n;
      V const r = rotate(q,p-o) + o;

      body.m_X(k) = r(0);
      body.m_Y(k) = r(1);
      body.m_Z(k) = r(2);
    }

    body.m_x = MT::convert(body.m_X, body.m_Y, body.m_Z);
  }


  template<typename MT>
  inline void bend(
                   Body<MT> & body
                   , typename MT::vector3_type const & line
                   , typename MT::vector3_type const & axis
                   , typename MT::real_type const & angle
                   )
  {
    typedef typename MT::real_type       T;
    typedef typename MT::vector3_type    V;
    typedef typename MT::value_traits    VT;
    typedef typename MT::quaternion_type Q;

    size_t const N = body.m_X.size();

    V const n = unit(line);

    V object_center  = V::make( body.m_X(0), body.m_Y(0), body.m_Z(0) );
    T w_start        = inner_prod(n, object_center );
    T w_end          = inner_prod(n, object_center );

    for(size_t k = 1u; k < N; ++k)
    {
      V const p = V::make( body.m_X(k), body.m_Y(k), body.m_Z(k) );
      T const w = inner_prod(n, p);

      w_start = w < w_start ? w : w_start;
      w_end   = w > w_end   ? w : w_end;

      object_center = object_center + p;
    }
    object_center = object_center / VT::numeric_cast( N );

    T const bend_radius = (w_end-w_start) / angle;
    V const bend_axis   = unit(axis);
    V const bend_center = object_center - bend_radius * cross(n, bend_axis);

    T const a_start =  -angle*VT::half();
    T const a_end   =  angle*VT::half();

    for(size_t k = 0u; k < N; ++k)
    {
      V const p = V::make( body.m_X(k), body.m_Y(k), body.m_Z(k) );
      T const w = inner_prod(n, p);
      T const t = (w-w_start) / (w_end-w_start);
      T const a = t*(a_end - a_start) + a_start;
      Q const q = Q::Ru(a, bend_axis);
      V const g = rotate(q,object_center - bend_center) + bend_center;
      V const h = object_center + w*n;
      V const s = rotate(q,p-h);

      body.m_X(k) = g(0)+s(0);
      body.m_Y(k) = g(1)+s(1);
      body.m_Z(k) = g(2)+s(2);
    }

    body.m_x = MT::convert(body.m_X, body.m_Y, body.m_Z);
  }

} // namespace hyper

// HYPER_MODIFIERS_H
#endif 
