#include <procedural_factory.h>

#include <eigenhelperall.h>

namespace procedural
{

template <typename T>
void make_heavy_sphere_light_sphere(content::API* engine, const EigenVector3<T>& position,
                                    const EigenQuaternion<T>& orientation, const T& sphere_radius, size_t const& scale,
                                    MaterialInfo<T> mat_info)
{

    T const bottom_stone_density = get_material_density_eigen<T>(mat_info, "Stone");
    T const top_stone_density = bottom_stone_density * scale;
    size_t const mid = get_material_id_eigen<T>(mat_info, "Stone");

    GeometryHandleEigen<T> sphere_handle = create_geometry_handle_sphere<T>(engine, sphere_radius);

    // light sphere, on ground
    const EigenVector3<T> T_b2m = sphere_handle.Tb2m();
    const EigenQuaternion<T> Q_b2m = sphere_handle.Qb2m();

    EigenVector3<T> T_m2l = EigenVector3<T>(0, sphere_radius, 0);
    const EigenQuaternion<T> Q_m2l = Rotateu(-std::numbers::pi_v<T> * 0.5f, EigenVector3<T>(1, 0, 0));

    const EigenVector3<T> T_l2w = (position);
    const EigenQuaternion<T> Q_l2w = (orientation);

    EigenVector3<T> T_b2w;
    EigenQuaternion<T> Q_b2w;

    compute_body_to_world_transform<T>(T_b2m, Q_b2m, T_m2l, Q_m2l, T_l2w, Q_l2w, T_b2w, Q_b2w);

    create_rigid_body<T>(engine, T_b2w, Q_b2w, sphere_handle, mid, bottom_stone_density);

    // heavy sphere on top of light sphere
    T_m2l = EigenVector3<T>(0, 3 * sphere_radius, 0);

    compute_body_to_world_transform<T>(T_b2m, Q_b2m, T_m2l, Q_m2l, T_l2w, Q_l2w, T_b2w, Q_b2w);

    create_rigid_body<T>(engine, T_b2w, Q_b2w, sphere_handle, mid, top_stone_density);
}

template void make_heavy_sphere_light_sphere<float>(content::API* engine, const EigenVector3<float>& position,
                                                    const EigenQuaternion<float>& orientation,
                                                    const float& sphere_radius, size_t const& scale,
                                                    MaterialInfo<float> mat_info);

} //namespace procedural
