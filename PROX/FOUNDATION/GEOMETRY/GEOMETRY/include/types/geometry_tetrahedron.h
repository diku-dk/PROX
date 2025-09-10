#ifndef GEOMETRY_TETRAHEDRON_H
#define GEOMETRY_TETRAHEDRON_H

#include <types/geometry_support_mapping.h>
#include <types/geometry_triangle.h>

#include <tiny_vector_functions.h>
#include <tiny_is_number.h>
#include <tiny_is_finite.h>

#include <cassert>
#include <cmath>     // needed for std::fabs

namespace geometry
{


template<typename T>
class TetrahedronEigen
    : public SupportMapping<T>
{
protected:
    EigenVector3<T> m_point[4];

public:

    const EigenVector3<T>& point(unsigned int const & idx ) const { return m_point[idx]; }
    EigenVector3<T>& point(unsigned int const & idx )       { return m_point[idx]; }
    const EigenVector3<T>& p(unsigned int const & idx )     const { return m_point[idx]; }
    EigenVector3<T>& p(unsigned int const & idx )           { return m_point[idx]; }

public:

    TetrahedronEigen()
    {
        m_point[0] = {0,0,0};
        m_point[1] = {1,0,0};
        m_point[2] = {0,1,0};
        m_point[3] = {0,0,1};
    }

    TetrahedronEigen(
        const EigenVector3<T>& p0
        , const EigenVector3<T>& p1
        , const EigenVector3<T>& p2
        , const EigenVector3<T>& p3
        )
    {
        m_point[0] = p0;
        m_point[1] = p1;
        m_point[2] = p2;
        m_point[3] = p3;
    }

    TetrahedronEigen( TetrahedronEigen const & tet)
    {
        *this = tet;
    }

    TetrahedronEigen & operator=(TetrahedronEigen const & tet)
    {
        if (this != &tet)
        {
            this->m_point[0] = tet.m_point[0];
            this->m_point[1] = tet.m_point[1];
            this->m_point[2] = tet.m_point[2];
            this->m_point[3] = tet.m_point[3];
        }
        return *this;
    }

public:
    EigenVector3<T> get_support_point(EigenVector3<T> v) const override
    {
        assert( is_number(v(0)) || !"get_support_point(): NAN encountered");
        assert( is_number(v(1)) || !"get_support_point(): NAN encountered");
        assert( is_number(v(2)) || !"get_support_point(): NAN encountered");
        assert( is_finite(v(0)) || !"get_support_point(): INF encountered");
        assert( is_finite(v(1)) || !"get_support_point(): INF encountered");
        assert( is_finite(v(2)) || !"get_support_point(): INF encountered");

        EigenVector3<T> p = this->m_point[0];

        assert( is_number(p(0)) || !"get_support_point(): NAN encountered");
        assert( is_number(p(1)) || !"get_support_point(): NAN encountered");
        assert( is_number(p(2)) || !"get_support_point(): NAN encountered");
        assert( is_finite(p(0)) || !"get_support_point(): INF encountered");
        assert( is_finite(p(1)) || !"get_support_point(): INF encountered");
        assert( is_finite(p(2)) || !"get_support_point(): INF encountered");

        T max_val = dot(p,v);

        for(size_t i =1u; i < 4u;)
        {
            const EigenVector3<T>& q = this->m_point[i++];

            assert( is_number(q(0)) || !"get_support_point(): NAN encountered");
            assert( is_number(q(1)) || !"get_support_point(): NAN encountered");
            assert( is_number(q(2)) || !"get_support_point(): NAN encountered");
            assert( is_finite(q(0)) || !"get_support_point(): INF encountered");
            assert( is_finite(q(1)) || !"get_support_point(): INF encountered");
            assert( is_finite(q(2)) || !"get_support_point(): INF encountered");

            T const val = dot(q,v);

            if(val > max_val)
            {
                max_val = val;
                p = q;
            }
        }
        assert( is_number(p(0)) || !"get_support_point(): NAN encountered");
        assert( is_number(p(1)) || !"get_support_point(): NAN encountered");
        assert( is_number(p(2)) || !"get_support_point(): NAN encountered");
        assert( is_finite(p(0)) || !"get_support_point(): INF encountered");
        assert( is_finite(p(1)) || !"get_support_point(): INF encountered");
        assert( is_finite(p(2)) || !"get_support_point(): INF encountered");

        return p;
    }

/*    EigenVector3<T> get_support_point(EigenVector3<T> dir) const override
    {
        auto bestPoint = m_point[0];
        auto maxVal = dot(dir, bestPoint);
        for (size_t i = 1; i < 4; ++i)
        {
            auto currPoint = m_point[i];
            auto currVal = dot(currPoint, dir);
            if (currVal > maxVal)
            {
                maxVal = currVal;
                bestPoint = currPoint;
            }
        }
        return bestPoint;
    }*/

    EigenVector3<T> get_center() const
    {
        const EigenVector3<T>& p0 = this->m_point[0];
        const EigenVector3<T>& p1 = this->m_point[1];
        const EigenVector3<T>& p2 = this->m_point[2];
        const EigenVector3<T>& p3 = this->m_point[3];

        return (p0 + p1 + p2 + p3)/ 4;
    }

    T get_scale() const
    {
        using std::max;
        using std::sqrt;

        const EigenVector3<T>& c  = this->get_center();
        const EigenVector3<T>& p0 = this->m_point[0];
        const EigenVector3<T>& p1 = this->m_point[1];
        const EigenVector3<T>& p2 = this->m_point[2];
        const EigenVector3<T>& p3 = this->m_point[3];

        T const dot0 = (p0-c).dot( p0-c);
        T const dot1 = (p1-c).dot( p1-c);
        T const dot2 = (p2-c).dot( p2-c);
        T const dot3 = (p3-c).dot( p3-c);

        T const max_val = max( dot0, max( dot1, max( dot2, dot3) ) );

        T const diameter = 2 * sqrt( max_val );

        assert(diameter > 0 || !"get_scale(): Internal error");
        assert(is_finite(diameter)   || !"get_scale(): Inf");
        assert(is_number(diameter)   || !"get_scale(): Nan");

        return diameter;
    }

};

template<typename T>
inline TetrahedronEigen<T> make_tetrahedron()
{
    return TetrahedronEigen<T>();
}

template<typename T>
inline TetrahedronEigen<T> make_tetrahedron(
    const EigenVector3<T>& p0
    , const EigenVector3<T>& p1
    , const EigenVector3<T>& p2
    , const EigenVector3<T>& p3
    )
{
    return TetrahedronEigen<T>(p0,p1,p2,p3);
}


template<typename T>
inline TetrahedronEigen<T> uniform_scale( const T& scale, TetrahedronEigen<T> const & A)
{
    const EigenVector3<T>& p0 = A.p(0u);
    const EigenVector3<T>& p1 = A.p(1u);
    const EigenVector3<T>& p2 = A.p(2u);
    const EigenVector3<T>& p3 = A.p(3u);
    const EigenVector3<T>& c  = A.get_center();

    return TetrahedronEigen<T>(
        (scale*(p0-c)) + c
        , (scale*(p1-c)) + c
        , (scale*(p2-c)) + c
        , (scale*(p3-c)) + c
        );
}

template<typename T>
inline T get_signed_volume(TetrahedronEigen<T> const & tet)
{

    const EigenVector3<T>& x1 = tet.point(0);
    const EigenVector3<T>& x2 = tet.point(1);
    const EigenVector3<T>& x3 = tet.point(2);
    const EigenVector3<T>& x4 = tet.point(3);

    T const vol6 = ( x4 - x1).dot( (x2 - x1).cross(x3 - x1) );

    assert( is_number(vol6) || !"get_signed_volume(): Nan encountered");
    assert( is_finite(vol6) || !"get_signed_volume(): Inf encountered");

    return vol6 / (6.0);
}

template<typename T>
inline T get_volume(TetrahedronEigen<T> const & tet)
{
    using std::fabs;

    return std::abs(get_signed_volume(tet));
}
}

namespace geometry
{

template <typename T>
inline Triangle<T> get_opposite_face(unsigned int const& idx,
                                     TetrahedronEigen<T> const& tet)
{
    assert(idx < 4u || !"get_opposite_face(): idx out of range");

    //
    // We assume that tetrahedron vertices are given with positive orientation
    //
    // That is vertex 3 is on the front side of triangle given by the ccw order
    // of vertices 0, 1, 2
    //
    //  face opposite vertex 0:  123
    //  face opposite vertex 1:  032
    //  face opposite vertex 2:  013
    //  face opposite vertex 3:  021

    if (idx == 0u) return Triangle<T>(tet.point(1), tet.point(2), tet.point(3));
    if (idx == 1u) return Triangle<T>(tet.point(0), tet.point(3), tet.point(2));
    if (idx == 2u) return Triangle<T>(tet.point(0), tet.point(1), tet.point(3));
    if (idx == 3u) return Triangle<T>(tet.point(0), tet.point(2), tet.point(1));

    return Triangle<T>();
}
} // namespace geometry

//GEOMETRY_TETRAHEDRON_H
#endif
