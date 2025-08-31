#include <convex_shapes.h>

#include <tiny_math_types.h>
#include <tiny_is_finite.h>

#include <cassert>
#include <cmath>    // needed for std::min and std::sqrt

namespace convex
{

  template<typename T>
  const T& Cylinder<T>::half_height() const { return this->m_half_height; }

  template<typename T>
  T& Cylinder<T>::half_height() { return this->m_half_height; }

  template<typename T>
  const T& Cylinder<T>::radius() const { return this->m_radius; }

  template<typename T>
  T& Cylinder<T>::radius() { return this->m_radius; }

  template<typename T>
  Cylinder<T>::Cylinder()
      : m_half_height( 1 )
      , m_radius( 1 )
  {}

  template<typename T>
  auto Cylinder<T>::get_support_point(EigenVector3<T> dir) const -> EigenVector3<T>
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

  template<typename T>
  T Cylinder<T>::get_scale() const
  {
    using std::min;

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

  template class Cylinder<float>;
  template class Cylinder<double>;

} // namespace convex
