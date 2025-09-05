#include <procedural_factory.h>

#include <tiny_math_types.h>

namespace procedural
{

template <typename T>
void make_wall(content::API* engine, const EigenVector3<T>& position, const EigenQuaternion<T>& orientation,
               const T& width, const T& height, const T& depth, size_t const& layers, size_t const& span,
               MaterialInfo<T> mat_info)
{

    T const stone_width = width / span;
    T const stone_height = height / layers;
    T const stone_density = get_material_density_eigen<T>(mat_info, "Stone");
    size_t const mid = get_material_id_eigen<T>(mat_info, "Stone");

    GeometryHandleEigen<T> stone_handle = create_geometry_handle_box<T>(engine, stone_width, stone_height, depth);

    for (size_t i = 0u; i < layers; ++i)
    {
        for (size_t j = 0u; j < span; ++j)
        {
            T const x = j * (stone_width) + ((i % 2) * 0.5f * stone_width);
            T const y = (i + 0.5f) * (stone_height);
            T const z = 0;

            const EigenVector3<T> T_b2m = stone_handle.Tb2m();
            const EigenQuaternion<T> Q_b2m = stone_handle.Qb2m();

            const EigenVector3<T> T_m2l = EigenVector3<T>(x, y, z);
            const EigenQuaternion<T> Q_m2l = EigenQuaternion<T>::Identity();

            const EigenVector3<T> T_l2w = position;
            const EigenQuaternion<T> Q_l2w = orientation;

            EigenVector3<T> T_b2w;
            EigenQuaternion<T> Q_b2w;

            compute_body_to_world_transform<T>(T_b2m, Q_b2m, T_m2l, Q_m2l, T_l2w, Q_l2w, T_b2w, Q_b2w);

            create_rigid_body<T>(engine, T_b2w, Q_b2w, stone_handle, mid, stone_density);
        }
    }
}

template void make_wall<float>(content::API* engine, const EigenVector3<float>& position,
                               const EigenQuaternion<float>& orientation, const float& width, const float& height,
                               const float& depth, size_t const& layers, size_t const& span,
                               MaterialInfo<float> mat_info);

} //namespace procedural
