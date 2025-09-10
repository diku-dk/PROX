#include <procedural_factory.h>

#include <eigenhelperall.h>

namespace procedural
{

template <typename T>
void make_ground(content::API* engine, const EigenVector3<T>& position, const EigenQuaternion<T>& orientation,
                 MaterialInfo<T> mat_info, const T& width, const T& height, const T& depth)
{
    size_t const mid = get_material_id_eigen<T>(mat_info, "Ground");

    GeometryHandleEigen<T> ground = create_geometry_handle_box<T>(engine, width, height, depth);

    const EigenVector3<T> Pw = rotate(orientation, EigenVector3<T>(0, -height * 0.5f, 0)) + position;
    const EigenQuaternion<T> Qw = orientation;

    create_rigid_body<T>(engine, Pw, Qw, ground, mid, 1, true, "Visualizer/ground");
}

template void make_ground<float>(content::API* engine, const EigenVector3<float>& position,
                                 const EigenQuaternion<float>& orientation, MaterialInfo<float> mat_info,
                                 const float& width, const float& height, const float& depth);

} //namespace procedural
