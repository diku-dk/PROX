#include <procedural_factory.h>

#include <tiny_math_types.h>

#include <cmath>

namespace procedural
{

template <typename T>
void make_arch(content::API* engine, const EigenVector3<T>& position, const EigenQuaternion<T>& orientation,
               const T& r_outer, const T& r_inner, const T& pillar_height, const T& stone_depth,
               size_t const& arch_slices, size_t const& pillar_segments, MaterialInfo<T> mat_info)
{
    using std::cos;
    using std::sin;

    T const delta_theta = std::numbers::pi_v<T> / arch_slices;
    T const pillar_stone_width = pillar_height / pillar_segments;
    T const pillar_stone_height = r_outer - r_inner;
    T const center_radius = (r_outer + r_inner) * 0.5f;
    T const stone_density = get_material_density_eigen<T>(mat_info, "Stone");
    size_t const mid = get_material_id_eigen<T>(mat_info, "Stone");

    std::vector<EigenVector3<T>> arch_vertices(8u);

    compute_arch_stone_vertices_eigen<T>(delta_theta, stone_depth, r_outer, r_inner, &arch_vertices[0]);

    GeometryHandleEigen<T> arch_stone = create_geometry_handle_cuboid<T>(engine, &arch_vertices[0]);
    GeometryHandleEigen<T> pillar_stone
        = create_geometry_handle_box<T>(engine, pillar_stone_width, pillar_stone_height, stone_depth);

    for (size_t i = 0; i < pillar_segments; ++i)
    {
        T const y = (0.5f + i) * pillar_stone_width;

        const EigenVector3<T> T_b2m = pillar_stone.Tb2m();
        const EigenQuaternion<T> Q_b2m = pillar_stone.Qb2m();

        const EigenVector3<T> T_left_m2l = EigenVector3<T>(-center_radius, y, 0);
        const EigenQuaternion<T> Q_left_m2l = Rotateu<T>(std::numbers::pi_v<T> * 0.5f, EigenVector3<T>({0, 0, 1}));

        const EigenVector3<T> T_l2w = position;
        const EigenQuaternion<T> Q_l2w = orientation;

        EigenQuaternion<T> Q_left_b2w;
        EigenVector3<T> T_left_b2w;

        compute_body_to_world_transform<T>(T_b2m, Q_b2m, T_left_m2l, Q_left_m2l, T_l2w, Q_l2w, T_left_b2w, Q_left_b2w);

        create_rigid_body<T>(engine, T_left_b2w, Q_left_b2w, pillar_stone, mid, stone_density);

        const EigenVector3<T> T_right_m2l = EigenVector3<T>(center_radius, y, 0);
        const EigenQuaternion<T> Q_right_m2l = Rotateu<T>(std::numbers::pi_v<T> * 0.5f, EigenVector3<T>({0, 0, 1}));

        EigenVector3<T> T_right_b2w;
        EigenQuaternion<T> Q_right_b2w;

        compute_body_to_world_transform<T>(T_b2m, Q_b2m, T_right_m2l, Q_right_m2l, T_l2w, Q_l2w, T_right_b2w,
                                           Q_right_b2w);

        create_rigid_body<T>(engine, T_right_b2w, Q_right_b2w, pillar_stone, mid, stone_density);
    }

    for (size_t i = 0; i < arch_slices; ++i)
    {
        T const theta = delta_theta * (i + 0.5f);
        T const x = center_radius * cos(theta);
        T const y = center_radius * sin(theta) + pillar_height;

        const EigenVector3<T> T_b2m = (arch_stone.Tb2m());
        const EigenQuaternion<T> Q_b2m = (arch_stone.Qb2m());

        const EigenVector3<T> T_m2l = EigenVector3<T>(x, y, 0);
        const EigenQuaternion<T> Q_m2l = Rotateu(theta - std::numbers::pi_v<T> * 0.5f, EigenVector3<T>({0, 0, 1}));
        const EigenVector3<T> T_l2w = position;
        const EigenQuaternion<T> Q_l2w = orientation;

        EigenVector3<T> T_b2w;
        EigenQuaternion<T> Q_b2w;

        compute_body_to_world_transform<T>(T_b2m, Q_b2m, T_m2l, Q_m2l, T_l2w, Q_l2w, T_b2w, Q_b2w);

        create_rigid_body<T>(engine, T_b2w, Q_b2w, arch_stone, mid, stone_density);
    }
}

template void make_arch(content::API* engine, const EigenVector3<float>& position,
                        const EigenQuaternion<float>& orientation, const float& r_outer, const float& r_inner,
                        const float& pillar_height, const float& stone_depth, size_t const& arch_slices,
                        size_t const& pillar_segments, MaterialInfo<float> mat_info);
} //namespace procedural
