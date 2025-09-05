#include <procedural_factory.h>

#include <tiny_math_types.h>

namespace procedural
{

template <typename T>
void make_stack(content::API* engine, const EigenVector3<T>& position, const EigenQuaternion<T>& orientation,
                const T& stone_dim, size_t const& layers, MaterialInfo<T> mat_info)
{
    T const stone_density = get_material_density_eigen<T>(mat_info, "Stone");
    size_t const mid = get_material_id_eigen<T>(mat_info, "Stone");

    GeometryHandleEigen<T> stone_handle = create_geometry_handle_box<T>(engine, stone_dim, stone_dim, stone_dim);

    for (size_t i = 0u; i < layers; ++i)
    {
        T y = (i + 0.5f) * stone_dim;

        const EigenVector3<T> T_b2m = stone_handle.Tb2m();
        const EigenQuaternion<T> Q_b2m = stone_handle.Qb2m();

        const EigenVector3<T> T_m2l = EigenVector3<T>(0, y, 0);
        const EigenQuaternion<T> Q_m2l = EigenQuaternion<T>::Identity();

        const EigenVector3<T> T_l2w = position;
        const EigenQuaternion<T> Q_l2w = orientation;

        EigenVector3<T> T_b2w;
        EigenQuaternion<T> Q_b2w;

        compute_body_to_world_transform<T>(T_b2m, Q_b2m, T_m2l, Q_m2l, T_l2w, Q_l2w, T_b2w, Q_b2w);

        create_rigid_body<T>(engine, T_b2w, Q_b2w, stone_handle, mid, stone_density);
    }
}

template void make_stack<float>(content::API* engine, const EigenVector3<float>& position,
                                const EigenQuaternion<float>& orientation, const float& stone_dim, size_t const& layers,
                                MaterialInfo<float> mat_info);

} //namespace procedural
