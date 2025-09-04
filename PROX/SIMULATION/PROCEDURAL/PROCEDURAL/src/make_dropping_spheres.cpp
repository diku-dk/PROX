#include <procedural_factory.h>

#include <tiny_math_types.h>

namespace procedural
{

template <typename T>
void make_dropping_spheres(content::API* engine, const EigenVector3<T>& position, const EigenQuaternion<T>& orientation,
                           const T& sphere_radius, size_t const& w, size_t const& h, size_t const& d,
                           MaterialInfo<T> mat_info)
{

    T const stone_density = get_material_density_eigen<T>(mat_info, "Stone");
    size_t const mid = get_material_id_eigen<T>(mat_info, "Stone");

    GeometryHandleEigen<T> sphere_handle = create_geometry_handle_sphere<T>(engine, sphere_radius);
    GeometryHandleEigen<T> plane_handle = create_geometry_handle_box<T>(engine, 2, (0.25f), 2);

    T x = 1;
    T y = 1;
    T z = 1;
    T const offset_x = sphere_radius * w;
    T const offset_z = sphere_radius * d;

    for (size_t i = 0; i < w; ++i)
    {
        x = 2 * i * sphere_radius - offset_x;
        for (size_t j = 0; j < h; ++j)
        {
            y = sphere_radius + 2 * j * sphere_radius + 4;

            for (size_t k = 0; k < d; ++k)
            {
                z = 2 * k * sphere_radius - offset_z;

                const EigenVector3<T> T_b2m = sphere_handle.Tb2m();
                const EigenQuaternion<T> Q_b2m = sphere_handle.Qb2m();

                const EigenVector3<T> T_m2l = EigenVector3<T>(x, y, z);
                const EigenQuaternion<T> Q_m2l = EigenQuaternion<T>::Identity();

                const EigenVector3<T> T_l2w = (position);
                const EigenQuaternion<T> Q_l2w = (orientation);

                EigenVector3<T> T_b2w;
                EigenQuaternion<T> Q_b2w;

                compute_body_to_world_transform<T>(T_b2m, Q_b2m, T_m2l, Q_m2l, T_l2w, Q_l2w, T_b2w, Q_b2w);

                create_rigid_body<T>(engine, T_b2w, Q_b2w, sphere_handle, mid, stone_density);
            }
        }
    }

    EigenVector3<T> T_b2m = plane_handle.Tb2m();
    EigenQuaternion<T> Q_b2m = plane_handle.Qb2m();

    EigenVector3<T> T_m2l = EigenVector3<T>((-1.25), (2.75), 0);
    EigenQuaternion<T> Q_m2l = Rotatez(-std::numbers::pi_v<T> * 0.25f);

    const EigenVector3<T> T_l2w = (position);
    const EigenQuaternion<T> Q_l2w = (orientation);

    EigenVector3<T> T_b2w;
    EigenQuaternion<T> Q_b2w;

    compute_body_to_world_transform<T>(T_b2m, Q_b2m, T_m2l, Q_m2l, T_l2w, Q_l2w, T_b2w, Q_b2w);

    create_rigid_body<T>(engine, T_b2w, Q_b2w, plane_handle, mid, stone_density, true);

    T_b2w = plane_handle.Tb2m();
    Q_b2w = plane_handle.Qb2m();

    T_m2l = EigenVector3<T>((1.25), (2.75), 0);
    Q_m2l = Rotatez(-std::numbers::pi_v<T> * 0.25f);

    compute_body_to_world_transform<T>(T_b2m, Q_b2m, T_m2l, Q_m2l, T_l2w, Q_l2w, T_b2w, Q_b2w);
    create_rigid_body<T>(engine, T_b2w, Q_b2w, plane_handle, mid, stone_density, true);
}

template void make_dropping_spheres<float>(content::API* engine, const EigenVector3<float>& position,
                                           const EigenQuaternion<float>& orientation, const float& sphere_radius,
                                           size_t const& w, size_t const& h, size_t const& d,
                                           MaterialInfo<float> mat_info);

} //namespace procedural
