#ifndef GEOMETRY_OVERLAP_OBB_CAPSULE
#define GEOMETRY_OVERLAP_OBB_CAPSULE

#include <types/geometry_obb.h>
#include <types/geometry_capsule.h>

#include <geometry_transform.h>

#include <eigenhelperall.h>

#include <cmath>
#include <vector>

namespace geometry
{
namespace detail
{

template <typename T>
inline void compute_obb_capsule_sat_axes(OBBEigen<T> const& obb,
                                         CapsuleEigen<T> const& capsule,
                                         std::vector<EigenVector3<T>>& axes)
{
    using std::fabs;
    using std::max;
    using std::min;

    const EigenMatrix3<T> R(obb.orientation());

    const EigenVector3<T> A0 = R.col(0);
    const EigenVector3<T> A1 = R.col(1);
    const EigenVector3<T> A2 = R.col(2);
    const EigenVector3<T> B = capsule.point1() - capsule.point0();

    assert(fabs(1 - dot(A0, A0)) < std::numeric_limits<T>::epsilon() * 10
           || !"compute_obb_capsule_sat_axes(): logic error");
    assert(fabs(1 - dot(A1, A1)) < std::numeric_limits<T>::epsilon() * 10
           || !"compute_obb_capsule_sat_axes(): logic error");
    assert(fabs(1 - dot(A2, A2)) < std::numeric_limits<T>::epsilon() * 10
           || !"compute_obb_capsule_sat_axes(): logic error");
    assert(fabs(dot(A0, A1)) < std::numeric_limits<T>::epsilon() * 10
           || !"compute_obb_capsule_sat_axes(): logic error");
    assert(fabs(dot(A0, A2)) < std::numeric_limits<T>::epsilon() * 10
           || !"compute_obb_capsule_sat_axes(): logic error");
    assert(fabs(dot(A1, A2)) < std::numeric_limits<T>::epsilon() * 10
           || !"compute_obb_capsule_sat_axes(): logic error");

    axes.resize(7u);

    axes[0] = A0;
    axes[1] = A1;
    axes[2] = A2;
    axes[3] = B;
    axes[4] = cross(A0, B);
    axes[5] = cross(A1, B);
    axes[6] = cross(A2, B);

    T const l3 = norm(axes[3]);
    T const l4 = norm(axes[4]);
    T const l5 = norm(axes[5]);
    T const l6 = norm(axes[6]);

    axes[3] = (l3 > std::numeric_limits<T>::epsilon() * 10) ? axes[3] / l3 : B;
    axes[4] = (l4 > std::numeric_limits<T>::epsilon() * 10) ? axes[4] / l4 : A0;
    axes[5] = (l5 > std::numeric_limits<T>::epsilon() * 10) ? axes[5] / l5 : A1;
    axes[6] = (l6 > std::numeric_limits<T>::epsilon() * 10) ? axes[6] / l6 : A2;

    assert(is_number(axes[0](0)) || !"compute_obb_capsule_sat_axes(): nan");
    assert(is_finite(axes[0](0)) || !"compute_obb_capsule_sat_axes(): inf");
    assert(is_number(axes[0](1)) || !"compute_obb_capsule_sat_axes(): nan");
    assert(is_finite(axes[0](1)) || !"compute_obb_capsule_sat_axes(): inf");
    assert(is_number(axes[0](2)) || !"compute_obb_capsule_sat_axes(): nan");
    assert(is_finite(axes[0](2)) || !"compute_obb_capsule_sat_axes(): inf");
    assert(fabs(1 - dot(axes[0], axes[0]))
               < std::numeric_limits<T>::epsilon() * 10
           || !"compute_obb_capsule_sat_axes(): logic error");

    assert(is_number(axes[1](0)) || !"compute_obb_capsule_sat_axes(): nan");
    assert(is_finite(axes[1](0)) || !"compute_obb_capsule_sat_axes(): inf");
    assert(is_number(axes[1](1)) || !"compute_obb_capsule_sat_axes(): nan");
    assert(is_finite(axes[1](1)) || !"compute_obb_capsule_sat_axes(): inf");
    assert(is_number(axes[1](2)) || !"compute_obb_capsule_sat_axes(): nan");
    assert(is_finite(axes[1](2)) || !"compute_obb_capsule_sat_axes(): inf");
    assert(fabs(1 - dot(axes[1], axes[1]))
               < std::numeric_limits<T>::epsilon() * 10
           || !"compute_obb_capsule_sat_axes(): logic error");

    assert(is_number(axes[2](0)) || !"compute_obb_capsule_sat_axes(): nan");
    assert(is_finite(axes[2](0)) || !"compute_obb_capsule_sat_axes(): inf");
    assert(is_number(axes[2](1)) || !"compute_obb_capsule_sat_axes(): nan");
    assert(is_finite(axes[2](1)) || !"compute_obb_capsule_sat_axes(): inf");
    assert(is_number(axes[2](2)) || !"compute_obb_capsule_sat_axes(): nan");
    assert(is_finite(axes[2](2)) || !"compute_obb_capsule_sat_axes(): inf");
    assert(fabs(1 - dot(axes[2], axes[2]))
               < std::numeric_limits<T>::epsilon() * 10
           || !"compute_obb_capsule_sat_axes(): logic error");

    assert(is_number(axes[3](0)) || !"compute_obb_capsule_sat_axes(): nan");
    assert(is_finite(axes[3](0)) || !"compute_obb_capsule_sat_axes(): inf");
    assert(is_number(axes[3](1)) || !"compute_obb_capsule_sat_axes(): nan");
    assert(is_finite(axes[3](1)) || !"compute_obb_capsule_sat_axes(): inf");
    assert(is_number(axes[3](2)) || !"compute_obb_capsule_sat_axes(): nan");
    assert(is_finite(axes[3](2)) || !"compute_obb_capsule_sat_axes(): inf");
    assert(fabs(1 - dot(axes[3], axes[3]))
               < std::numeric_limits<T>::epsilon() * 10
           || !"compute_obb_capsule_sat_axes(): logic error");

    assert(is_number(axes[4](0)) || !"compute_obb_capsule_sat_axes(): nan");
    assert(is_finite(axes[4](0)) || !"compute_obb_capsule_sat_axes(): inf");
    assert(is_number(axes[4](1)) || !"compute_obb_capsule_sat_axes(): nan");
    assert(is_finite(axes[4](1)) || !"compute_obb_capsule_sat_axes(): inf");
    assert(is_number(axes[4](2)) || !"compute_obb_capsule_sat_axes(): nan");
    assert(is_finite(axes[4](2)) || !"compute_obb_capsule_sat_axes(): inf");
    assert(fabs(1 - dot(axes[4], axes[4]))
               < std::numeric_limits<T>::epsilon() * 10
           || !"compute_obb_capsule_sat_axes(): logic error");

    assert(is_number(axes[5](0)) || !"compute_obb_capsule_sat_axes(): nan");
    assert(is_finite(axes[5](0)) || !"compute_obb_capsule_sat_axes(): inf");
    assert(is_number(axes[5](1)) || !"compute_obb_capsule_sat_axes(): nan");
    assert(is_finite(axes[5](1)) || !"compute_obb_capsule_sat_axes(): inf");
    assert(is_number(axes[5](2)) || !"compute_obb_capsule_sat_axes(): nan");
    assert(is_finite(axes[5](2)) || !"compute_obb_capsule_sat_axes(): inf");
    assert(fabs(1 - dot(axes[5], axes[5]))
               < std::numeric_limits<T>::epsilon() * 10
           || !"compute_obb_capsule_sat_axes(): logic error");

    assert(is_number(axes[6](0)) || !"compute_obb_capsule_sat_axes(): nan");
    assert(is_finite(axes[6](0)) || !"compute_obb_capsule_sat_axes(): inf");
    assert(is_number(axes[6](1)) || !"compute_obb_capsule_sat_axes(): nan");
    assert(is_finite(axes[6](1)) || !"compute_obb_capsule_sat_axes(): inf");
    assert(is_number(axes[6](2)) || !"compute_obb_capsule_sat_axes(): nan");
    assert(is_finite(axes[6](2)) || !"compute_obb_capsule_sat_axes(): inf");
    assert(fabs(1 - dot(axes[6], axes[6]))
               < std::numeric_limits<T>::epsilon() * 10
           || !"compute_obb_capsule_sat_axes(): logic error");
}

}// end namespace detail

  /**
   *
   */
template <typename T>
inline bool overlap_obb_capsule(OBBEigen<T> const& obb,
                                CapsuleEigen<T> const& capsule)
{
    using std::fabs;
    using std::max;
    using std::min;

    std::vector<EigenVector3<T>> axes;

    detail::compute_obb_capsule_sat_axes(obb, capsule, axes);

    std::vector<EigenVector3<T>> a(8u, EigenVector3<T>(0, 0, 0));
    a[0] = transform_from_obb(get_local_corner(0, obb), obb);
    a[1] = transform_from_obb(get_local_corner(1, obb), obb);
    a[2] = transform_from_obb(get_local_corner(2, obb), obb);
    a[3] = transform_from_obb(get_local_corner(3, obb), obb);
    a[4] = transform_from_obb(get_local_corner(4, obb), obb);
    a[5] = transform_from_obb(get_local_corner(5, obb), obb);
    a[6] = transform_from_obb(get_local_corner(6, obb), obb);
    a[7] = transform_from_obb(get_local_corner(7, obb), obb);

    std::vector<T> a_min(7u, std::numeric_limits<T>::max());
    std::vector<T> b_min(7u, std::numeric_limits<T>::max());
    std::vector<T> a_max(7u, std::numeric_limits<T>::lowest());
    std::vector<T> b_max(7u, std::numeric_limits<T>::lowest());

    T min_overlap = std::numeric_limits<T>::lowest();

    for (size_t i = 0u; i < 7u; ++i)
    {

        for (typename std::vector<EigenVector3<T>>::const_iterator p_a
             = a.begin();
             p_a != a.end(); ++p_a)
        {
            T const d = dot((*p_a), axes[i]);
            a_min[i] = min(a_min[i], d);
            a_max[i] = max(a_max[i], d);
        }

        {
            T const d0 = dot(capsule.point0(), axes[i]);
            T const d1 = dot(capsule.point1(), axes[i]);

            b_min[i] = min(b_min[i], d0 - capsule.radius());
            b_max[i] = max(b_max[i], d0 + capsule.radius());

            b_min[i] = min(b_min[i], d1 - capsule.radius());
            b_max[i] = max(b_max[i], d1 + capsule.radius());
        }

        assert(is_number(a_min[i]) || !"overlap_obb_capsule(); nan");
        assert(is_finite(a_min[i]) || !"overlap_obb_capsule(): inf");
        assert(is_number(a_max[i]) || !"overlap_obb_capsule(); nan");
        assert(is_finite(a_max[i]) || !"overlap_obb_capsule(): inf");

        assert(is_number(b_min[i]) || !"overlap_obb_capsule(); nan");
        assert(is_finite(b_min[i]) || !"overlap_obb_capsule(): inf");
        assert(is_number(b_max[i]) || !"overlap_obb_capsule(); nan");
        assert(is_finite(b_max[i]) || !"overlap_obb_capsule(): inf");

        if (a_max[i] < b_min[i]) return false;

        if (b_max[i] < a_min[i]) return false;

        if (a_min[i] <= b_min[i] && b_min[i] <= a_max[i])
        {
            T const overlap = b_min[i] - a_max[i];
            if (overlap > min_overlap)
            {
                min_overlap = overlap;
//          n = axes[i];
            }
        }
        if (b_min[i] <= a_min[i] && a_min[i] <= b_max[i])
        {
            T const overlap = a_min[i] - b_max[i];
            if (overlap > min_overlap)
            {
                min_overlap = overlap;
//          n = -axes[i];
            }
        }
        assert(is_number(min_overlap) || !"overlap_obb_capsule(); nan");
        assert(is_finite(min_overlap) || !"overlap_obb_capsule(): inf");
    }

    return (min_overlap <= 0);
}

} //namespace geometry

  // GEOMETRY_OVERLAP_OBB_CAPSULE
#endif
