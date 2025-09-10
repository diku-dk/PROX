#include <procedural.h>
#include <procedural_factory.h>

#include <eigenhelperall.h>

namespace procedural
{

template <typename T>
void make_pantheon(content::API* engine, const EigenVector3<T>& position, const EigenQuaternion<T>& orientation,
                   const T& r_outer, const T& r_inner, const T& height, size_t const& slices, size_t const& segments,
                   MaterialInfo<T> mat_info)
{
    const EigenQuaternion<T> Qw = orientation;
    EigenVector3<T> pw = rotate(Qw, position);

    make_tower<T>(engine, pw, Qw, r_outer, r_inner, height, slices * 2, segments, mat_info, false);

    pw = pw + rotate(Qw, EigenVector3<T>(0, 0, height));

    make_dome<T>(engine, pw, Qw, r_outer, r_inner, slices, segments, mat_info);
}

template void make_pantheon<float>(content::API* engine, const EigenVector3<float>& position,
                                   const EigenQuaternion<float>& orientation, const float& r_outer,
                                   const float& r_inner, const float& height, size_t const& slices,
                                   size_t const& segments, MaterialInfo<float> mat_info);

} //namespace procedural
