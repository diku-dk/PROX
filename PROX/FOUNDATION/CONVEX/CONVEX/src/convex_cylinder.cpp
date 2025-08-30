#include <convex_shapes.h>

#include <tiny_math_types.h>
#include <tiny_is_finite.h>

#include <cassert>
#include <cmath>    // needed for std::min and std::sqrt

namespace convex
{

  template<typename M>
  typename M::real_type const & Cylinder<M>::half_height() const { return this->m_half_height; }

  template<typename M>
  typename M::real_type & Cylinder<M>::half_height() { return this->m_half_height; }

  template<typename M>
  typename M::real_type const & Cylinder<M>::radius() const { return this->m_radius; }

  template<typename M>
  typename M::real_type & Cylinder<M>::radius() { return this->m_radius; }

  template<typename M>
  Cylinder<M>::Cylinder()
  : m_half_height( 1 )
  , m_radius( 1 )
  {}

  template<typename M>
  auto Cylinder<M>::get_support_point(EigenVector3<T> dir) const -> EigenVector3<T>
  {
      auto height = dir.z() > 0 ? m_half_height : -m_half_height;
      auto sigmaNorm = std::sqrt(dir.x()*dir.x() + dir.y()*dir.y());
      if (sigmaNorm > 0)
      {
          return
          {
              m_radius * dir.x() / sigmaNorm,
              m_radius * dir.y() / sigmaNorm,
              height
          };
      }
      return {0, 0, height};
  }

  template<typename M>
  typename M::real_type Cylinder<M>::get_scale() const
  {
    using std::min;

    typedef typename M::real_type     T;
    typedef typename M::value_traits VT;

    assert( is_number(this->m_half_height)    || !"NAN encountered");
    assert( is_finite(this->m_half_height)    || !"INF encountered");
    assert( this->m_half_height >= 0 || !"Negative half height");
    assert( is_number(this->m_radius)         || !"NAN encountered");
    assert( is_finite(this->m_radius)         || !"INF encountered");
        assert( this->m_radius >= 0      || !"Negative radius");

    T const d = 2 * ((this->m_radius > 0) ? this->m_radius : std::numeric_limits<T>::max());
    T const h = 2 * ((this->m_half_height > 0) ? this->m_half_height : std::numeric_limits<T>::max());

    return min(h, d);
  }

  using Mf = tiny::MathTypes<float>;
  using Md = tiny::MathTypes<double>;

  template class Cylinder<Mf>;
  template class Cylinder<Md>;

} // namespace convex
