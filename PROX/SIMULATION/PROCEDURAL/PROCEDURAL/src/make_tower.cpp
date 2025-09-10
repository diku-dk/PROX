#include <procedural_factory.h>

#include <eigenhelperall.h>

namespace procedural
{

template <typename T>
void make_tower(content::API* engine, const EigenVector3<T>& position, const EigenQuaternion<T>& orientation,
                const T& r_outer, const T& r_inner, const T& height, size_t const& slices, size_t const& segments,
                MaterialInfo<T> mat_info, bool const& use_cubes)
{

    size_t const mid = get_material_id_eigen<T>(mat_info, "Stone");

    T const stone_depth = height / segments;
    T const delta_theta = 2 * std::numbers::pi_v<T> / slices;
    T const center_radius = (r_outer + r_inner) * 0.5f;
    T const stone_density = get_material_density_eigen<T>(mat_info, "Stone");

    std::vector<EigenVector3<T>> vertices(8u);
    compute_arch_stone_vertices_eigen<T>(delta_theta, stone_depth, r_outer, r_inner, &vertices[0]);

    //T const center_height    = (r_outer + r_inner)*0.5f;
    T const max_height = r_outer * cos(delta_theta * 0.5f);
    T const min_height = r_inner * cos(delta_theta * 0.5f);
    T const half_width_inner = r_inner * sin(delta_theta * 0.5f);
    //T const half_width_outer = r_outer * sin( delta_theta*0.5f );
    T const half_depth = stone_depth * 0.5f;

    GeometryHandleEigen<T> stone_handle;

    if (use_cubes)
    {
        stone_handle
            = create_geometry_handle_box<T>(engine, 1.9 * half_width_inner, max_height - min_height, 2 * half_depth);
    }
    else { stone_handle = (create_geometry_handle_cuboid<T>(engine, &vertices[0])); }

    for (size_t i = 0u; i < segments; ++i)
    {
        for (size_t j = 0u; j < slices; ++j)
        {
            T const theta = j * delta_theta + (i % 2) * (delta_theta * 0.5f);
            T const x = center_radius * cos(theta);
            T const y = center_radius * sin(theta);
            T const z = (i + 0.5f) * stone_depth;

            const EigenVector3<T> T_b2m = stone_handle.Tb2m();
            const EigenQuaternion<T> Q_b2m = stone_handle.Qb2m();

            const EigenVector3<T> T_m2l = EigenVector3<T>(x, y, z);
            const EigenQuaternion<T> Q_m2l = Rotateu<T>(theta - std::numbers::pi_v<T> * 0.5, EigenVector3<T>(0, 0, 1));

            const EigenVector3<T> T_l2w = position;
            const EigenQuaternion<T> Q_l2w = orientation;

            EigenVector3<T> T_b2w;
            EigenQuaternion<T> Q_b2w;

            compute_body_to_world_transform<T>(T_b2m, Q_b2m, T_m2l, Q_m2l, T_l2w, Q_l2w, T_b2w, Q_b2w);

            create_rigid_body<T>(engine, T_b2w, Q_b2w, stone_handle, mid, stone_density);
        }
    }
}

template void make_tower<float>(content::API* engine, const EigenVector3<float>& position,
                                const EigenQuaternion<float>& orientation, const float& r_outer, const float& r_inner,
                                const float& height, size_t const& slices, size_t const& segments,
                                MaterialInfo<float> mat_info, bool const& use_cubes);

} //namespace procedural
