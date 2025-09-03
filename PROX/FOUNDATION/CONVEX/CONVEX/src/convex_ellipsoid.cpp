#include <convex_shapes.h>

#include <tiny_math_types.h>
#include <tiny_is_finite.h>

#include <cassert>
#include <cmath>     // needed for std::min and std::sqrt

namespace convex
{

template <typename T> const EigenVector3<T>& Ellipsoid<T>::scale() const { return this->m_scale; }

template <typename T> void Ellipsoid<T>::setScale(const EigenVector3<T>& vec) { this->m_scale = vec; }

template <typename T> const EigenVector3<T>& Ellipsoid<T>::scale() { return this->m_scale; }

template <typename T>
Ellipsoid<T>::Ellipsoid()
    : m_scale(1)
{
}

template <typename T> auto Ellipsoid<T>::get_support_point(EigenVector3<T> dir) const -> EigenVector3<T>
{
    T const& vx = dir(0);
    T const& vy = dir(1);
    T const& vz = dir(2);

    assert(is_number(vx) || !"NAN encountered");
    assert(is_number(vy) || !"NAN encountered");
    assert(is_number(vz) || !"NAN encountered");
    assert(is_finite(vx) || !"INF encountered");
    assert(is_finite(vy) || !"INF encountered");
    assert(is_finite(vz) || !"INF encountered");

    T const& sx = this->m_scale(0);
    T const& sy = this->m_scale(1);
    T const& sz = this->m_scale(2);

    assert(sx || !"NAN encountered");
    assert(sx || !"INF encountered");
    assert(sy || !"NAN encountered");
    assert(sy || !"INF encountered");
    assert(sz || !"NAN encountered");
    assert(sz || !"INF encountered");
    assert(sx >= 0 || !"Negative scale encountered");
    assert(sy >= 0 || !"Negative scale encountered");
    assert(sz >= 0 || !"Negative scale encountered");

      /*
    // An ellipsoid, E, is simply a scaled unit ball, B, and a scale is a linear
    // transformation, T. We can write it in a general way as
    //
    //  E = T(B)
    //
    // That means we can create a support function, S, of an ellipsoid from
    // that of a unit sphere ball by
    //
    //  S_E(v) = S_{T(B)}(v)
    //
    // Further we know that for any affine transformation, T(v) = R v + t, we have
    //
    //   S_{T(B)}(v)  =   T(  S_B( R^T v )  )
    //
    // In our particular case R = R^T = D, where D = diag(s_0,s_1,s_3) and t=0. Here
    // the s_i's are the axes scales repsectively. Putting it all together we have
    //
    //   S_{E}(v)  =   D(  S_B( D v )  )
    //
    // This is the formula implemented by this functor.
    */

    T const vv = vx * vx + vy * vy + vz * vz;

    assert(is_number(vv) || !"NAN encountered");
    assert(is_finite(vv) || !"INF encountered");

    T px = 0;
    T py = 0;
    T pz = 0;

    if (vv > 0)
    {
        T const wx = vx * sz;
        T const wy = vy * sy;
        T const wz = vz * sz;
        T const ww = wx * wx + wy * wy + wz * wz;

        assert(is_number(ww) || !"NAN encountered");
        assert(is_finite(ww) || !"INF encountered");

        T const tmp = 1 / sqrt(ww);

        assert(is_number(tmp) || !"NAN encountered");
        assert(is_finite(tmp) || !"INF encountered");

        px = wx * tmp * sx;
        py = wy * tmp * sy;
        pz = wz * tmp * sz;
    }
    else
    {
        px = sx;
        py = 0;
        pz = 0;
    }

    assert(is_number(px) || !"NAN encountered");
    assert(is_number(py) || !"NAN encountered");
    assert(is_number(pz) || !"NAN encountered");
    assert(is_finite(px) || !"INF encountered");
    assert(is_finite(py) || !"INF encountered");
    assert(is_finite(pz) || !"INF encountered");

    return {px, py, pz};
}

template <typename T> T Ellipsoid<T>::get_scale() const
{
    using std::min;

    T const& sx = this->m_scale(0);
    T const& sy = this->m_scale(1);
    T const& sz = this->m_scale(2);

    assert(is_number(sx) || !"NAN encountered");
    assert(is_finite(sx) || !"INF encountered");
    assert(is_number(sy) || !"NAN encountered");
    assert(is_finite(sy) || !"INF encountered");
    assert(is_number(sz) || !"NAN encountered");
    assert(is_finite(sz) || !"INF encountered");
    assert(sx >= 0 || !"Negative scale encountered");
    assert(sy >= 0 || !"Negative scale encountered");
    assert(sz >= 0 || !"Negative scale encountered");

    T const w = 2 * ((sx > 0) ? sx : std::numeric_limits<T>::max());
    T const h = 2 * ((sy > 0) ? sy : std::numeric_limits<T>::max());
    T const d = 2 * ((sz > 0) ? sz : std::numeric_limits<T>::max());

    return min(w, min(h, d));
}
template class Ellipsoid<float>;
template class Ellipsoid<double>;

} // namespace convex
