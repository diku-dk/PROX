#include <convex_shapes.h>

#include <types/geometry_sphere.h>

#include <tiny_math_types.h>
#include <tiny_is_finite.h>

#include <cassert>
#include <cmath>   // needed for std::min and std::sqrt

namespace convex
{

  template<typename M>
  typename M::real_type const & Capsule<M>::half_height() const { return this->m_half_height; }

  template<typename M>
  typename M::real_type & Capsule<M>::half_height() { return this->m_half_height; }

  template<typename M>
  typename M::real_type const & Capsule<M>::radius() const { return this->m_radius; }

  template<typename M>
  typename M::real_type & Capsule<M>::radius() { return this->m_radius; }

  template<typename M>
  Capsule<M>::Capsule()
  : m_half_height( 1 )
  , m_radius( 1 )
  {}

  template<typename M>
  auto Capsule<M>::get_support_point(EigenVector3<T> dir) const -> EigenVector3<T>
  {
      using std::sqrt;

            assert( is_number(this->m_half_height)    || !"NAN encountered");
      assert( is_finite(this->m_half_height)    || !"INF encountered");
      assert( this->m_half_height >= 0 || !"Negative half height");
      assert( is_number(this->m_radius)         || !"NAN encountered");
      assert( is_finite(this->m_radius)         || !"INF encountered");
      assert( this->m_radius >= 0      || !"Negative radius");

      geometry::Sphere<typename M::vector3_type> S;
      S.radius() = this->m_radius;

      // Get the support point of the sphere
      EigenVector3<T> p = S.get_support_point(dir);

      // Cut the sphere into two halves and displace them along the z-axis.
      if( dir(2) > 0 )
          p(2) += this->m_half_height;
      else if( dir(2) < 0 )
          p(2) -= this->m_half_height;


      return p;
  }

  template<typename M>
  typename M::real_type Capsule<M>::get_scale() const
  {
    using std::min;

    typedef typename M::real_type    T;
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

  template class Capsule<Mf>;
  template class Capsule<Md>;

} // namespace convex
