#include <procedural_factory.h>

#include <tiny_math_types.h>

#include <cmath>

namespace procedural
{

template <typename T>
void make_box(content::API* engine, const EigenVector3<T>& position, const EigenQuaternion<T>& orientation,
              const T& width, const T& height, const T& depth, MaterialInfo<T> mat_info, bool const fixed)
{
    using std::cos;
    using std::sin;

    T const stone_density = get_material_density_eigen<T>(mat_info, "Stone");
    size_t const mid = get_material_id_eigen<T>(mat_info, "Stone");

    GeometryHandleEigen<T> box_handle = create_geometry_handle_box<T>(engine, width, height, depth);

    const EigenVector3<T> T_b2m = box_handle.Tb2m();
    const EigenQuaternion<T> Q_b2m = box_handle.Qb2m();

    const EigenVector3<T> T_m2l = EigenVector3<T>(0, 0, 0);
    const EigenQuaternion<T> Q_m2l = EigenQuaternion<T>::Identity();

    const EigenVector3<T> T_l2w = position;
    const EigenQuaternion<T> Q_l2w = orientation;

    EigenVector3<T> T_b2w;
    EigenQuaternion<T> Q_b2w;

    compute_body_to_world_transform<T>(T_b2m, Q_b2m, T_m2l, Q_m2l, T_l2w, Q_l2w, T_b2w, Q_b2w);

    create_rigid_body<T>(engine, T_b2w, Q_b2w, box_handle, mid, stone_density, fixed);
}

template void make_box<float>(content::API* engine, const EigenVector3<float>& position,
                              const EigenQuaternion<float>& orientation, const float& width, const float& height,
                              const float& depth, MaterialInfo<float> mat_info, bool const fixed);

} //namespace procedural
