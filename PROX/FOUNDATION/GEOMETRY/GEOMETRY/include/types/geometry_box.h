#ifndef GEOMETRY_BOX_H
#define GEOMETRY_BOX_H

#include <types/geometry_support_mapping.h>

#include <cassert>

namespace geometry
{

  template<typename V>
  class Box
  : public geometry::SupportMapping<V>
  {
  public:

    typedef typename V::real_type       T;
    typedef typename V::value_traits    VT;

  protected:

    V m_half_extent;

  public:

    V const & half_extent() const   { return this->m_half_extent; }
    V       & half_extent()         { return this->m_half_extent; }

  public:

    Box()
    : m_half_extent( V::make(1,1, 1) )
    {}

    Box(V const &  half_extent )
    : m_half_extent( half_extent )
    {}

    Box(Box const & box)
    {
      *this = box;
    }

    Box & operator=(Box const & box)
    {
      if( this != &box)
      {
        this->m_half_extent = box.m_half_extent;
      }
      return *this;
    }

  public:
    EigenVector3<T> get_support_point(EigenVector3<T> v) const override
    {
        EigenVector3<T> result;
        for (size_t i = 0; i < 3; ++i)
        {
            result[i] = v[i] > 0 ? m_half_extent(i) : -m_half_extent(i);
        }
        return result;
    }

    T get_scale() const
    {
      using std::min;

      T const & hx = this->m_half_extent(0);
      T const & hy = this->m_half_extent(1);
      T const & hz = this->m_half_extent(2);

      assert( is_number(hx)    || !"NAN encountered");
      assert( is_finite(hx)    || !"INF encountered");
      assert( hx >= 0 || !"Negative half extent encountered");
      assert( is_number(hy)    || !"NAN encountered");
      assert( is_finite(hy)    || !"INF encountered");
      assert( hy >= 0 || !"Negative half extent encountered");
      assert( is_number(hz)    || !"NAN encountered");
      assert( is_finite(hz)    || !"INF encountered");
      assert( hz >= 0 || !"Negative half extent encountered");

      T const w = 2 * ((hx > 0) ? hx : std::numeric_limits<T>::max());
      T const h = 2 * ((hy > 0) ? hy : std::numeric_limits<T>::max());
      T const d = 2 * ((hz > 0) ? hz : std::numeric_limits<T>::max());

      return min( w, min(h, d) );
    }

  };

  template<typename V>
  inline Box<V> make_box(V const &  half_extent )
  {
    return Box<V>( half_extent);
  }

  /**
   * Returns the coordiantes of a specifed corner point.
   */
  template<typename V>
  inline V get_corner(size_t const & idx, Box<V> const & box)
  {
    assert(idx < 8u || !"get_local_corner(): Idx was out of range");

    typedef typename V::real_type    T;

    T const x = ((idx   ) & 0x0001u) ? box.half_extent()(0) : -box.half_extent()(0);
    T const y = ((idx>>1) & 0x0001u) ? box.half_extent()(1) : -box.half_extent()(1);
    T const z = ((idx>>2) & 0x0001u) ? box.half_extent()(2) : -box.half_extent()(2);

    return V::make(x,y,z);
  }

} //namespace geometry

// GEOMETRY_BOX_H
#endif
