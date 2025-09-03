#include <convex_shapes.h>

#include <types/geometry_sphere.h>

#include <tiny_math_types.h>
#include <tiny_is_finite.h>

#include <cassert>
#include <cmath>   // needed for std::min and std::sqrt

namespace convex
{

template <typename T> const T& Capsule<T>::half_height() const { return this->m_half_height; }

template <typename T> T& Capsule<T>::half_height() { return this->m_half_height; }

template <typename T> const T& Capsule<T>::radius() const { return this->m_radius; }

template <typename T> T& Capsule<T>::radius() { return this->m_radius; }

template <typename T>
Capsule<T>::Capsule()
    : m_half_height(1)
    , m_radius(1)
{
}

template <typename T> auto Capsule<T>::get_support_point(EigenVector3<T> dir) const -> EigenVector3<T>
{
    using std::sqrt;

    assert(is_number(this->m_half_height) || !"NAN encountered");
    assert(is_finite(this->m_half_height) || !"INF encountered");
    assert(this->m_half_height >= 0 || !"Negative half height");
    assert(is_number(this->m_radius) || !"NAN encountered");
    assert(is_finite(this->m_radius) || !"INF encountered");
    assert(this->m_radius >= 0 || !"Negative radius");

    geometry::Sphere<T> S;
    S.radius() = this->m_radius;

      // Get the support point of the sphere
    EigenVector3<T> p = S.get_support_point(dir);

      // Cut the sphere into two halves and displace them along the z-axis.
    if (dir(2) > 0)
        p(2) += this->m_half_height;
    else if (dir(2) < 0)
        p(2) -= this->m_half_height;

    return p;
}

template <typename T> T Capsule<T>::get_scale() const
{
    using std::min;

    assert(is_number(this->m_half_height) || !"NAN encountered");
    assert(is_finite(this->m_half_height) || !"INF encountered");
    assert(this->m_half_height >= 0 || !"Negative half height");
    assert(is_number(this->m_radius) || !"NAN encountered");
    assert(is_finite(this->m_radius) || !"INF encountered");
    assert(this->m_radius >= 0 || !"Negative radius");

    T const d = 2 * ((this->m_radius > 0) ? this->m_radius : std::numeric_limits<T>::max());
    T const h = 2 * ((this->m_half_height > 0) ? this->m_half_height : std::numeric_limits<T>::max());

    return min(h, d);
}
template class Capsule<float>;
template class Capsule<double>;

} // namespace convex
