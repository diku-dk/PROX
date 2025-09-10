#include <convex_shapes.h>

#include <eigenhelperall.h>
#include <cassert>
#include <cmath>    // needed for std::min

namespace convex
{

template <typename T> void ConvexHull<T>::add_point(const EigenVector3<T>& p)
{
    assert(is_number(p(0)) || !"nan encountered");
    assert(is_number(p(1)) || !"nan encountered");
    assert(is_number(p(2)) || !"nan encountered");
    assert(is_finite(p(0)) || !"infinite encountered");
    assert(is_finite(p(1)) || !"infinite encountered");
    assert(is_finite(p(2)) || !"infinite encountered");

    this->m_points.push_back(p);
}

template <typename T> const EigenVector3<T>& ConvexHull<T>::get_point(size_t const& idx) const
{
    assert(idx < this->size() || !"get_point(): illegal index");
    return this->m_points[idx];
}

template <typename M> size_t ConvexHull<M>::size() const { return this->m_points.size(); }

template <typename M> void ConvexHull<M>::clear() { this->m_points.clear(); }

template <typename T>
ConvexHull<T>::ConvexHull()
    : m_points()
{
}

template <typename T> auto ConvexHull<T>::get_support_point(EigenVector3<T> v) const -> EigenVector3<T>
{

    assert(is_number(v(0)) || !"NAN encountered");
    assert(is_number(v(1)) || !"NAN encountered");
    assert(is_number(v(2)) || !"NAN encountered");
    assert(is_finite(v(0)) || !"INF encountered");
    assert(is_finite(v(1)) || !"INF encountered");
    assert(is_finite(v(2)) || !"INF encountered");

    size_t const N = this->m_points.size();
    assert(N > 0u || !"empty hull");

    EigenVector3<T> p = (this->m_points[0]);

    assert(is_number(p(0)) || !"NAN encountered");
    assert(is_number(p(1)) || !"NAN encountered");
    assert(is_number(p(2)) || !"NAN encountered");
    assert(is_finite(p(0)) || !"INF encountered");
    assert(is_finite(p(1)) || !"INF encountered");
    assert(is_finite(p(2)) || !"INF encountered");

    T max_val = p.dot(v);

    for (size_t i = 1u; i < N;)
    {
        EigenVector3<T> q = (this->m_points[i++]);

        assert(is_number(q(0)) || !"NAN encountered");
        assert(is_number(q(1)) || !"NAN encountered");
        assert(is_number(q(2)) || !"NAN encountered");
        assert(is_finite(q(0)) || !"INF encountered");
        assert(is_finite(q(1)) || !"INF encountered");
        assert(is_finite(q(2)) || !"INF encountered");

        T const val = q.dot(v);

        if (val > max_val)
        {
            max_val = val;
            p = q;
        }
    }
    assert(is_number(p(0)) || !"NAN encountered");
    assert(is_number(p(1)) || !"NAN encountered");
    assert(is_number(p(2)) || !"NAN encountered");
    assert(is_finite(p(0)) || !"INF encountered");
    assert(is_finite(p(1)) || !"INF encountered");
    assert(is_finite(p(2)) || !"INF encountered");

    return p;
}

template <typename T> T ConvexHull<T>::get_scale() const
{
    using std::min;

    size_t const N = this->m_points.size();

    assert(N > 0u || !"empty hull");

    T max_x = this->m_points[0](0);
    T max_y = this->m_points[0](1);
    T max_z = this->m_points[0](2);

    assert(is_number(max_x) || !"NAN encountered");
    assert(is_number(max_y) || !"NAN encountered");
    assert(is_number(max_z) || !"NAN encountered");
    assert(is_finite(max_x) || !"INF encountered");
    assert(is_finite(max_y) || !"INF encountered");
    assert(is_finite(max_z) || !"INF encountered");

    T min_x = max_x;
    T min_y = max_y;
    T min_z = max_z;

    for (size_t i = 1u; i < N;)
    {
        const EigenVector3<T>& q = this->m_points[i++];

        T const& x = q(0);
        T const& y = q(1);
        T const& z = q(2);

        assert(is_number(x) || !"NAN encountered");
        assert(is_number(y) || !"NAN encountered");
        assert(is_number(z) || !"NAN encountered");
        assert(is_finite(x) || !"INF encountered");
        assert(is_finite(y) || !"INF encountered");
        assert(is_finite(z) || !"INF encountered");

        max_x = (max_x > x) ? max_x : x;
        max_y = (max_y > y) ? max_y : y;
        max_z = (max_z > z) ? max_z : z;

        min_x = (min_x < x) ? min_x : x;
        min_y = (min_y < y) ? min_y : y;
        min_z = (min_z < z) ? min_z : z;
    }

    assert(is_number(max_x) || !"NAN encountered");
    assert(is_number(max_y) || !"NAN encountered");
    assert(is_number(max_z) || !"NAN encountered");
    assert(is_finite(max_x) || !"INF encountered");
    assert(is_finite(max_y) || !"INF encountered");
    assert(is_finite(max_z) || !"INF encountered");

    assert(is_number(min_x) || !"NAN encountered");
    assert(is_number(min_y) || !"NAN encountered");
    assert(is_number(min_z) || !"NAN encountered");
    assert(is_finite(min_x) || !"INF encountered");
    assert(is_finite(min_y) || !"INF encountered");
    assert(is_finite(min_z) || !"INF encountered");

    assert(max_x > min_x || !"Degenerate hull");
    assert(max_y > min_y || !"Degenerate hull");
    assert(max_z > min_z || !"Degenerate hull");

    T const w = max_x - min_x;
    T const h = max_y - min_y;
    T const d = max_z - min_z;

    return min(w, min(h, d));
}

template class ConvexHull<float>;
template class ConvexHull<double>;

} // namespace convex
